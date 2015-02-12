.. highlight:: c++

.. default-domain:: cpp

.. _documentation-image:

=====
Image
=====

Theia provides a basic :class:`Image\<T\>` class to use for SfM and SLAM
applications. The class is fairly lightweight, and is largely a wrapper for the
`CImg <http://cimg.sourceforge.net/>`_ but allows for more general use. It is
templated on the data type for the pixels (e.g. float, uchar, int) and is the
standard image class that used for feature detection within Theia.

.. class:: Image<T>

  Images can be stored either as RGB or grayscale, and the number of channels
  used to describe an image can be determined with `:func:Channels`. We have
  `typdef Image<float> FloatImage` for convenience. Images can be easily read in
  from a file upon contruction or with the :func:`Read` command.

  .. code-block:: c++

    FloatImage my_rgb_image(img_filename);

    // Convert to grayscale.
    FloatImage my_grayscale_image = my_rgb_image.AsGrayscaleImage();

    // ... do some stuff ... //

    // Write the image back out.
    my_rgb_image.Write(out_filename);
    my_grayscale_image.Write(out_filename);

  Once an image is loaded, pixel values can be accessed with ``()``
  operators. Pixels are referenced in x, y, c order (c is the color channel).

  .. code-block:: c++

    // Load RGB and scale images.
    FloatImage my_img("test_img.jpg");

    // Get the middle pixel location.
    int middle_y = my_img.Rows()/2;
    int middle_x = my_img.Cols()/2;

    // Grab the middle pixel.
    const float middle_r_pixel = my_img(middle_x, middle_y, 0);
    const float middle_g_pixel = my_img(middle_x, middle_y, 1);
    const float middle_b_pixel = my_img(middle_x, middle_y, 2);

    // Output the RGB Pixel value.
    LOG(INFO) << "red = " << middle_r_pixel.red
              << " green = " << middle_g_pixel.green
              << " blue = " << middle_b_pixel.blue;

    // Convert to grayscale.
    my_img.ConvertToGrayscaleImage();

    // Get the grayscale pixel value.
    const float middle_gray_pixel = my_img(middle_x, middle_y);

    // Output the grayscale pixel value.
    LOG(INFO) << "gray = " << middle_gray_pixel;

  We also read in the EXIF information to determine the focal length. Currently,
  we only make the focal length publicly accesible but we save all EXIF data and
  more functionality can be easily added if needed.

  .. code-block:: c++

    FloatImage my_img("test_img.jpg");
    const double focal_length_in_pixels;
    if (my_img.FocalLengthPixels(&focal_length_pixels)) {
      LOG(INFO) << "Focal length from EXIF is: " << focal_length_pixels;
    } else {
      LOG(INFO) << "Could not extract the focal length from EXIF data.";
    }

We have also implemented some useful member functions of the :class:`Image` class. For a full list of functions, `theia/image/image.h`

.. function:: int Image\<T\>::Rows() const
.. function:: int Image\<T\>::Cols() const
.. function:: int Image\<T\>::Channels() const
.. function:: T* Image\<T\>::Data()
.. function:: const T* Image\<T\>::Data() const
.. function:: void Image\<T\>::Read(const std::string& filename)
.. function:: void Image\<T\>::Write(const std::string& filename)
.. function:: void Image\<T\>::ConvertToGrayscaleImage()
.. function:: void Image\<T\>::ConvertToRGBImage()
.. function:: Image<T> Image\<T\>::AsGrayscaleImage() const
.. function:: Image<T> Image\<T\>::AsRGBImage() const
.. function:: Image\<T\> Image\<T\>::Integrate() const
.. function:: void Image\<T\>::Resize(int new_rows, int new_cols)
.. function:: void Image\<T\>::Resize(double scale)
.. function:: void Image\<T\>::HalfSample(Image\<T\>* out_image) const
.. function:: void Image\<T\>::TwoThirdsSample(Image\<T\>* out_image) const
