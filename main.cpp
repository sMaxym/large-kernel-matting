#include <iostream>
#include <string>
#include <fstream>

#include "include/cfmatting/ImageMatrix.h"
#include "include/cfmatting/CGMattingSolver.h"
#include "include/io/imageIO.h"

int main(int argc, char* argv[])
{
	std::string file_path = "image.jpg";
	FILE *in_file;
	if ((in_file = fopen(file_path.c_str(), "rb")) == nullptr) {
		std::cout << "Cannot open file " << file_path << std::endl;
		return -1;
	}
	ImageMatrix im = read_jpeg(in_file);
	fclose(in_file);

	file_path = "trimap.jpg";
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

	CGMattingSolver solver(im, trimap);
	auto alpha = solver.alphaMatting(10);

	std::ofstream out_fs("output_new.txt", std::fstream::out);

	for (size_t row = 0; row < im.height(); ++row)
	{
		for (size_t col = 0; col < im.width(); ++col)
		{
			out_fs << alpha[im.width() * row + col] << " ";
		}
		out_fs << std::endl;
	}

	out_fs.close();

	return 0;
}