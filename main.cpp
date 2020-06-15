#include <iostream>
#include <string>

#include "include/matting/ImageMatrix.h"
#include "include/matting/CGMattingSolver.h"
#include "include/io/imageIO.h"

int main(int argc, char* argv[])
{
	std::string file_path = "image_apple.jpg";
	FILE *in_file;
	if ((in_file = fopen(file_path.c_str(), "rb")) == nullptr) {
		std::cout << "Cannot open file " << file_path << std::endl;
		return -1;
	}
	ImageMatrix im = read_jpeg(in_file);
	fclose(in_file);

	file_path = "trimap_apple.jpg";
	if ((in_file = fopen(file_path.c_str(), "rb")) == nullptr) {
		std::cout << "Cannot open file " << file_path << std::endl;
		return -1;
	}
	ImageMatrix trimap_buffer = read_jpeg(in_file);
	fclose(in_file);

	ImageMatrix trimap(trimap_buffer.width(), trimap_buffer.height(), 1);
	for (size_t row = 0; row < trimap_buffer.height(); ++row)
	{
		for (size_t col = 0; col < trimap_buffer.width(); ++col)
		{
			trimap.setAt(row, col, 0, trimap_buffer.getAt(row, col, 0));
		}
	}

	im.normalize();
	trimap.normalize();

	CGMattingSolver solver(im, trimap, 20);
	solver.setRegParameter(0.001);
	auto alpha = solver.alphaMatting(20);

	ImageMatrix alpha_image(im.width(), im.height(), 3);
	for (size_t row = 0; row < im.height(); ++row)
	{
		for (size_t col = 0; col < im.width(); ++col)
		{
			alpha_image.setAt(row, col, 0, alpha[im.width() * row + col]);
			alpha_image.setAt(row, col, 1, alpha[im.width() * row + col]);
			alpha_image.setAt(row, col, 2, alpha[im.width() * row + col]);
		}
	}

	alpha_image.toImageFormat();

	if ((in_file = fopen("output.jpg", "wb")) == nullptr) {
		std::cout << "Cannot create file " << file_path << std::endl;
		return -1;
	}

	write_jpeg(in_file, alpha_image, 100);

	fclose(in_file);

	return 0;
}