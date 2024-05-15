#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <random>
std::default_random_engine generator;
std::normal_distribution<double> N01(0., 1.);

void sliced(double* src, double* tgt, int Npix) {
	std::vector<std::pair<double, int>> projection_src(Npix);
	std::vector<double> projection_tgt(Npix);

	for (int i = 0; i < 50; i++) {
		double xLine = N01(generator);
		double yLine = N01(generator);
		double zLine = N01(generator);
		double norm = sqrt(xLine * xLine + yLine * yLine + zLine * zLine);
		xLine /= norm;
		yLine /= norm;
		zLine /= norm;
		for (int j = 0; j < Npix; j++) {
			projection_src[j] = std::pair<double, int>(src[j * 3] * xLine + src[j * 3 + 1] * yLine + src[j * 3 + 2] * zLine, j);
			projection_tgt[j] = tgt[j * 3] * xLine + tgt[j * 3 + 1] * yLine + tgt[j * 3 + 2] * zLine;
		}
		std::sort(projection_src.begin(), projection_src.end());
		std::sort(projection_tgt.begin(), projection_tgt.end());
		for (int j = 0; j < Npix; j++) {
			src[projection_src[j].second * 3] += (projection_tgt[j] - projection_src[j].first) * xLine;
			src[projection_src[j].second * 3 + 1] += (projection_tgt[j] - projection_src[j].first) * yLine;
			src[projection_src[j].second * 3 + 2] += (projection_tgt[j] - projection_src[j].first) * zLine;
		}
	}
}

int main() {

	int W, H, C;
	int W2, H2, C2;
	
	//stbi_set_flip_vertically_on_load(true);
	unsigned char *image = stbi_load("8733654151_b9422bb2ec_k.jpg",
                                 &W,
                                 &H,
                                 &C,
                                 STBI_rgb);
	unsigned char *target = stbi_load("redim.jpg",
                                 &W2,
                                 &H2,
                                 &C2,
                                 STBI_rgb);
	std::vector<double> image_double(W*H*3);
	std::vector<double> target_double(W2*H2*3);

	for (int i=0; i<W*H*3; i++)
		image_double[i] = image[i];
	for (int i=0; i<W2*H2*3; i++)
		target_double[i] = target[i];

	sliced(&image_double[0], &target_double[0], W*H);

	std::vector<unsigned char> image_result(W*H * 3, 0);
	for (int i = 0; i < W2*H2*3; i++) {
		image_result[i] = (unsigned char) std::max(0., std::min(255., image_double[i]));
	}
	stbi_write_png("image.png", W, H, 3, &image_result[0], 0);

	return 0;
}