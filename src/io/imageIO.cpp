#include "./../../include/io/imageIO.h"

// TODO add exception handlers
ImageMatrix read_jpeg(FILE *jpeg_file)
{
	const size_t KBuffer_height = 1;
	size_t row_size;
	JSAMPROW buffer[KBuffer_height];
	struct jpeg_decompress_struct dcomp{};
	struct jpeg_error_mgr jerr{};
	dcomp.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&dcomp);

	jpeg_stdio_src(&dcomp, jpeg_file);
	jpeg_read_header(&dcomp, true);
	row_size = dcomp.image_width * dcomp.num_components;
	buffer[0] = new JSAMPLE[row_size];
	ImageMatrix image(dcomp.image_width, dcomp.image_height, dcomp.num_components);

	jpeg_start_decompress(&dcomp);
	while (dcomp.output_scanline < dcomp.image_height)
	{
		jpeg_read_scanlines(&dcomp, buffer, KBuffer_height);
		for (size_t col = 0; col < dcomp.image_width; ++col)
		{
			for (uint8_t channel = 0, val; channel < dcomp.num_components; ++channel)
			{
				val = buffer[0][dcomp.num_components * col + channel];
				image.setAt(dcomp.output_scanline - 1, col, channel, val);
			}
		}
	}

	jpeg_finish_decompress(&dcomp);
	jpeg_destroy_decompress(&dcomp);
	delete[] buffer[0];
	return image;
}

void write_jpeg(FILE *jpeg_file, const ImageMatrix& image, const size_t quality)
{
	struct jpeg_compress_struct compressor{};
	struct jpeg_error_mgr jerr{};
	compressor.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&compressor);
	jpeg_stdio_dest(&compressor, jpeg_file);

	compressor.image_width = image.width();
	compressor.image_height = image.height();
	compressor.input_components = image.numComponents();
	compressor.in_color_space = JCS_RGB;
	jpeg_set_defaults(&compressor);
	jpeg_set_quality(&compressor, quality, true);

	jpeg_start_compress(&compressor, true);
	JSAMPROW buffer[1];
	size_t row_size;
	row_size = image.width() * image.numComponents();

	buffer[0] = new JSAMPLE[row_size];
	while (compressor.next_scanline < compressor.image_height)
	{
		for (size_t col = 0; col < image.width(); ++col)
		{
			for (uint8_t channel = 0, val; channel < image.numComponents(); ++channel)
			{
				val = image.getAt(compressor.next_scanline, col, channel);
				val = std::lround(val);
				buffer[0][image.numComponents() * col + channel] = val;
			}
		}
		jpeg_write_scanlines(&compressor, buffer, 1);
	}
	jpeg_finish_compress(&compressor);
	jpeg_destroy_compress(&compressor);
	delete[] buffer[0];
}

ImageMatrix img_read(const std::string& path)
{
	FILE *data_file;
	if ((data_file = fopen(path.c_str(), "rb")) == nullptr)
	{
		throw std::runtime_error("cannot find jpeg");
	}

	ImageMatrix img;
	try
	{
		img = read_jpeg(data_file);
	} catch (const std::exception& e)
	{
		fclose(data_file);
		throw;
	}
	fclose(data_file);
	return img;
}

void img_write(const std::string& path, const ImageMatrix& img, int quality)
{
	FILE *data_file;
	if ((data_file = fopen(path.c_str(), "wb")) == nullptr)
	{
		throw std::runtime_error("cannot find or create file");
	}
	try
	{
		write_jpeg(data_file, img, quality);
	} catch (const std::exception& e)
	{
		fclose(data_file);
		throw;
	}
	fclose(data_file);
}