#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

#include "include/matting/ImageMatrix.h"
#include "include/matting/CGMattingSolver.h"
#include "include/io/imageIO.h"
#include "include/timer/timer.h"

void sigmoid_norming(ImageMatrix& image, double slope);

int main(int argc, char* argv[])
{
	const size_t KKernel_radius = 20,
				 KIterations = 20,
				 KOutput_quality = 100;
	const double KReg_param = 0.01;

	if (argc != 4)
	{
		std::cout << "wrong number of arguments" << std::endl;
		return 0;
	}

	std::string image_path = argv[1],
				trimap_path = argv[2],
				output_path= argv[3];

	ImageMatrix im, trimap;
	try
	{
		std::cout << "Image loading..." << std::endl;
		im = img_read(image_path);
		trimap = img_read(trimap_path).grayscale();
	} catch (const std::exception& e)
	{
		std::cout << e.what() << std::endl;
		return 0;
	}

	im.normalize();
	trimap.normalize();

	std::cout << "Solver init..." << std::endl;
	auto start_time = get_current_time_fenced();
	Eigen::VectorXd alpha;
	CGMattingSolver solver(im, trimap, KKernel_radius);
	solver.setRegParameter(KReg_param);
	try
	{
		std::cout << "Matting..." << std::endl;
		alpha = solver.alphaMatting(KIterations);
	} catch (const std::exception& e)
	{
		std::cout << e.what() << std::endl;
		return 0;
	}
	auto cur_time = get_current_time_fenced();
	std::cout << "Matting time total: " << to_us(cur_time - start_time)
	          << "ms" << std::endl;

	ImageMatrix alpha_image(alpha, im.width(), im.height());
	alpha_image.normalize();

	sigmoid_norming(alpha_image, 10);

	alpha_image.expandColorspace();
	alpha_image.expandColorspace();

	alpha_image.toImageFormat();
	try
	{
		std::cout << "Alpha mate saving..." << std::endl;
		img_write(output_path, alpha_image, KOutput_quality);
	} catch (const std::exception& e)
	{
		std::cout << e.what() << std::endl;
		return 0;
	}

	return 0;
}

void sigmoid_norming(ImageMatrix& image, const double slope)
{
	for (size_t row = 0; row < image.height(); ++row)
	{
		for (size_t col = 0; col < image.width(); ++col)
		{
			double value = image.getAt(row, col, 0);
			value = 1 / (1 + std::exp(-slope * (value - 0.5)));
			image.setAt(row, col, 0, value);
		}
	}
}

