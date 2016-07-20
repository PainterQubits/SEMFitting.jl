module SEMFitting
using Images
import DataFrames
using DataFrames: DataFrame, sort!
using Clustering

export approxscale!
export sizerange
export colrow
export groundplane

"""
```
approxscale!(img::AbstractImage, scale, scalerow=1758)
```

Given an SEM image `img` with scale bar of physical length `scale`,
returns physical length per pixel and sets the `spacedirections` property
of the image. Assumes that the scale bar is alone in pixel row `scalerow` of the
image, which defaults to 1758 (good for Painter group SEM).

Currently uses the far extrema of the scale bar to determine length.
"""
function approxscale!(img::AbstractImage, scale, scalerow=1758)
    gray = convert(Image{Images.Gray}, img)
    tf = convert(Array{Bool}, data(gray[:,scalerow]))
    last = tf[1]
    flag1 = true
    local firstedge
    local lastedge
    for (i,x) in enumerate(tf[2:end])
        if x == true && last == false && flag1
            firstedge = i
            flag1 = false
        elseif x == false && last == true
            lastedge = i-1
        end
        last = x
    end
    distpp = scale/(lastedge-firstedge)
    img["spacedirections"] = [distpp,0,0,distpp]
    distpp
end

"""
```
sizerange(img::AbstractImage, scale, min, max; thresh=0.5, tone=true)
```

Does grayscale conversion and thresholding of image `img`, with `thresh` for
thresholding and `tone` for optional inversion (default to true).

Next, extract components of the binary image that are within a specified size
range `min` to `max` where these are units of area in terms of the approximate
physical scale of the image. The `scale` parameter is the length of the scale bar.

Returns `im, df` where `im` is a black and white image only showing the
objects in the specified size range, and `df` is a `DataFrame` with the sizes,
image indices, and centroids of the found objects.
"""
function sizerange(img::AbstractImage, scale, min, max; thresh=0.5, tone=true)
    distpp = approxscale!(img,scale)
    a = convert(Image{Images.Gray}, img)
    a[a.<thresh] = 0.0
    a[a.>=thresh] = 1.0
    tone && (a = 1.0 - a)
    b = convert(Array{Bool}, Images.data(a))
    arr = label_components(b)
    df = DataFrames.DataFrame(len=component_lengths(arr).*(distpp^2), ind=component_indices(arr),
        cen=component_centroids(arr))
    df[:cen] = map(x->(x[1]*distpp, x[2]*distpp), df[:cen])
    sort!(df, cols=[:len])
    similar(arr, Bool)
    inds = reduce(vcat, df[min .<= df[:len] .<= max, :ind])
    c = similar(arr, Bool)
    c[:] = false
    c[inds] = true
    im = grayim(convert(Array{Float64}, c))
    im["spacedirections"] = img["spacedirections"]
    im, df[min .<= df[:len] .<= max, :]
end

"""
```
colrow(df::DataFrame, cols, rows)
```

Uses some clustering to determine rows and columns of found objects in `df`.
You should tell it the number of columns and rows for the correct clustering.
"""
function colrow(df::DataFrame, cols, rows)
    xlen = Array{Float64}(size(df)[1])
    ylen = Array{Float64}(size(df)[1])
    for (i,v) in enumerate(df[:cen])
        x,y = v
        xlen[i] = x
        ylen[i] = y
    end
    xcen = kmeans(Array{Float64}(transpose(xlen)), cols).centers
    ycen = kmeans(Array{Float64}(transpose(ylen)), rows).centers
    sort!(reshape(transpose(xcen), (length(xcen),))),
        sort!(reshape(transpose(ycen), (length(ycen),)))
    # transpose(xcen), transpose(ycen)
end

# img has correct scaling set
"""
```
groundplane(img::AbstractImage, cols_or_rows, dim, spacing, noiselen)
```

Given some image without rotation with rectangular holes in a ground plane...

- `img` - a binary image (expects grayscale).
- `cols_or_rows`- typically one of the results from [`colrow`](@ref).
- `dim` - 1 or 2 (should be 1 for col, 2 for row)
- `spacing` - The design spacing between hole centers in the ground plane.
- `noiselen` - Contiguous lengths of 1 or 0 shorter than this are ignored in
length determination.
"""
function groundplane(img::AbstractImage, cols_or_rows, dim, spacing, noiselen)
    sd = img["spacedirections"]
    distpp = sd[1]
    a = convert(Image{Images.Gray}, img)
    b = convert(Array{Bool}, Images.data(a))

    pixlengths = Array{Array{Int,1},1}(length(cols_or_rows))
    # pixlengthsrow = Array{Array{Int,1},1}(length(rows))

    for (i,x) in enumerate(cols_or_rows)
        ctr = 1
        flag = false
        pix = Array{Int}(0)
        pixlengths[i] = Array{Int}(0)
        cut = (dim == 1 ? b[Int(round(x/distpp)),:] : b[:, Int(round(x/distpp))])
        last = cut[1]
        for z in cut[2:end]
            if z != last
                push!(pix, ctr)
                ctr = 0
            end
            ctr += 1
            last = z
        end
        ctr != 1 && push!(pix, ctr)

        minpix = round(noiselen/distpp)
        runninglength = 0
        for v in pix
            runninglength += v
            if runninglength > minpix
                push!(pixlengths[i], runninglength)
                runninglength = 0
            end
        end
    end

    l = 0.0
    for arr in pixlengths
        deleteat!(arr, 1)
        pop!(arr)
        mod(length(arr), 2) == 1 && pop!(arr)
        l += reduce(+, arr) / (length(arr) / 2)
    end
    l /= length(pixlengths)
    # Now l is the average length in pixels
    distpp = spacing/l

    display(pixlengths)
    physlengths = reduce(hcat, [map(x->x*distpp, arr) for arr in pixlengths])
    # distpp, physlengths
end

end # module
