#ifndef CFIMAGEMATTING_IMAGEIO_H
#define CFIMAGEMATTING_IMAGEIO_H

#include <cstdio>
#include <jpeglib.h>
#include "./../cfmatting/ImageMatrix.h"

ImageMatrix read_jpeg(FILE *jpeg_file);
void write_jpeg(FILE *jpeg_file, const ImageMatrix& image, size_t quality=5);

#endif
