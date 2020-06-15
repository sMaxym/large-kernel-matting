#ifndef CFIMAGEMATTING_IMAGEIO_H
#define CFIMAGEMATTING_IMAGEIO_H

#include <cstdio>
#include <iostream>
#include <jpeglib.h>
#include "./../matting/ImageMatrix.h"

ImageMatrix read_jpeg(FILE *jpeg_file);
void write_jpeg(FILE *jpeg_file, const ImageMatrix& image, size_t quality=5);

ImageMatrix img_read(const std::string& path);
void img_write(const std::string& path, const ImageMatrix& img, int quality);

#endif
