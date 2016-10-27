.. highlight:: c++

.. default-domain:: cpp

.. _documentation-image:

=====
Image
=====

Theia provides a basic :class:`FloatImage` class to use for SfM and SLAM
applications. The class is fairly lightweight, and is largely a wrapper for the
`OpenImageIO <http://openimageio.org//>`_ but allows for more general use. It is
the standard image class that used for feature detection within Theia.

.. class:: FloatImage

  Images can be stored either as RGB or grayscale, and the number of channels
  used to describe an image can be determined with `:func:Channels`. Images can
  be easily read in from a file upon contruction or with the :func:`Read`
  command.

  .. code-block:: c++

    FloatImage my_rgb_image(img_filename);

    // Convert to grayscale.
    FloatImage my_grayscale_image = my_rgb_image.AsGrayscaleImage();

    // ... do some stuff ... //

    // Write the image back out.
    my_rgb_image.Write(out_filename);
    my_grayscale_image.Write(out_filename);

  Once an image is loaded, pixel values can be accessed through several accessor
  functions. To avoid ambiguity between matrix and image notation, we provide
  explicit methods to get image pixel values based on either the xy coordinate
  or the row and column.

  .. code-block:: c++

    // Load RGB and scale images.
    FloatImage my_img("test_img.jpg");

    // Get the middle pixel location.
    int middle_x = my_img.Width()/2;
    int middle_y = my_img.Height()/2;

    // Grab the middle pixel.
    const Eigen::Vector3f middle_xy_pixel = my_img.GetXY(middle_x, middle_y);

    // Get the middle pixel location.
    int middle_row = my_img.Rows()/2;
    int middle_col = my_img.Cols()/2;

    // Grab the middle pixel.
    const Eigen::Vector3f middle_row_pixel = my_img.GetRowCol(middle_row, middle_col);

    // The rgb values should be the same for both pixels.
    CHECK_EQ(middle_xy_pixel(0), middle_row_pixel(0));
    CHECK_EQ(middle_xy_pixel(1), middle_row_pixel(1));
    CHECK_EQ(middle_xy_pixel(2), middle_row_pixel(2));

  Similar operations apply for grayscale images and for obtaining the value of
  specific color channels at a pixel. The :func:`BilinearInterpolate` function
  allows for interpolation at non-discrete pixel locations. See the image.h file
  for more details.

We have also implemented some useful member functions of the :class:`FloatImage` class. For a full list of functions, `theia/image/image.h`

.. function:: int FloatImage::Rows() const
.. function:: int FloatImage::Cols() const
.. function:: int FloatImage::Channels() const
.. function:: float* FloatImage::Data()
.. function:: const float* FloatImage::Data() const
.. function:: void FloatImage::Read(const std::string& filename)
.. function:: void FloatImage::Write(const std::string& filename)
.. function:: void FloatImage::ConvertToGrayscaleImage()
.. function:: void FloatImage::ConvertToRGBImage()
.. function:: FloatImage FloatImage::AsGrayscaleImage() const
.. function:: FloatImage FloatImage::AsRGBImage() const
.. function:: FloatImage FloatImage::Integrate() const
.. function:: FloatImage FloatImage::ComputeGradient() const
.. function:: void FloatImage::Resize(int new_width, int new_height)
.. function:: void FloatImage::ResizeRowsCols(int new_rows, int new_cols)
.. function:: void FloatImage::Resize(double scale)
