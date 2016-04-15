
<a id='SEMFitting.jl-1'></a>

# SEMFitting.jl


A package for fitting images from the SEM in the Painter Group cleanroom.


<a id='Installation-1'></a>

## Installation


  * `Pkg.clone("https://github.com/ajkeller34/SEMFitting.jl.git")`


<a id='Usage-1'></a>

## Usage


Forthcoming


<a id='API-1'></a>

## API

<a id='SEMFitting.approxscale!' href='#SEMFitting.approxscale!'>#</a>
**`SEMFitting.approxscale!`** &mdash; *Function*.

---


`approxscale!(img::AbstractImage, scale, scalerow=1758)`

Given an SEM image `img` with scale bar of physical length `scale`, returns physical length per pixel and sets the `spacedirections` property of the image. Assumes that the scale bar is alone in pixel row `scalerow` of the image, which defaults to 1758 (good for Painter group SEM).

Currently uses the far extrema of the scale bar to determine length.

<a id='SEMFitting.colrow' href='#SEMFitting.colrow'>#</a>
**`SEMFitting.colrow`** &mdash; *Function*.

---


`colrow(df::DataFrame, cols, rows)`

Uses some clustering to determine rows and columns of found objects in `df`. You should tell it the number of columns and rows for the correct clustering.

<a id='SEMFitting.sizerange' href='#SEMFitting.sizerange'>#</a>
**`SEMFitting.sizerange`** &mdash; *Function*.

---


`sizerange(img::AbstractImage, scale, min, max; thresh=0.5, tone=true)`

Does grayscale conversion and thresholding of image `img`, with `thresh` for thresholding and `tone` for optional inversion (default to true).

Next, extract components of the binary image that are within a specified size range `min` to `max` where these are units of area in terms of the approximate physical scale of the image. The `scale` parameter is the length of the scale bar.

Returns `im, df` where `im` is a black and white image only showing the objects in the specified size range, and `df` is a `DataFrame` with the sizes, image indices, and centroids of the found objects.

<a id='SEMFitting.groundplane' href='#SEMFitting.groundplane'>#</a>
**`SEMFitting.groundplane`** &mdash; *Function*.

---


`groundplane(img::AbstractImage, cols_or_rows, dim, spacing, noiselen, minlen, maxlen)`

Given some image without rotation with rectangular holes in a ground plane...

  * `img` - a binary image (expects grayscale).
  * `cols_or_rows`- typically one of the results from [`colrow`](index.md#SEMFitting.colrow).
  * `dim` - 1 or 2 (should be 1 for col, 2 for row)
  * `spacing` - The design spacing between hole centers in the ground plane.
  * `noiselen` - Contiguous lengths of 1 or 0 shorter than this are ignored in length determination.
  * `minlen` - Reject lengths shorter than this from the final results.
  * `maxlen` - Reject lengths longer than this from the final results.

