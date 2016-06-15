/*
 Compiling:
 A simple console command 'make' will excute the compilation instructions
 found in the Makefile bundled with this code. It will result in an executable
 named filters.
 
 Running:
 The executable filters takes commands of the form:
 ./filters problem_ID input_File ouput_File additional_Arguments
 
 Some examples:
 
 ./filters avgGray input.png out.png
 
 This runs the 'averageGrayscale' function defined in this file and described
 in the project doc Part 2, problem 1 on the input. The output is saved as out.png.
 
 ----
 
 ./filters gaussFilter input.jpg output.png 0.5
 
 This runs the Gaussian filter function (gaussianFilter) with sigma = 0.5 on input.jpg.
 This is problem 2b from Part 4 in the project documentation.
 
 Test git commit to see the changes
 
 */
 
//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <math.h>
#include <algorithm>
 
//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;
 
//Basic image operations
CImg<double> averageGrayscale(CImg<double> input);
CImg<double> simpleBW(CImg<double> input);
CImg<double> advancedBW(CImg<double> input);
 
//Adding noise
CImg<double> uniformNoise(CImg<double> input);
CImg<double> gaussianNoise(CImg<double> input, double sigma);
CImg<double> saltAndPepperNoise(CImg<double> input);
 
//Filtering
CImg<double> filter(CImg<double> input, CImg<double> filter);
CImg<double> meanFilter(CImg<double> input, int filterSize);
CImg<double> gaussianFilter(CImg<double> input, double sigma);
CImg<double> medianFilter(CImg<double> input, int size);

//Analysis
double meanSquaredError(CImg<double> original, CImg<double> filtered);

//Optimization
CImg<double> separableKernelFilter(CImg<double> input, CImg<double>  filterX, CImg<double>  filterY);
CImg<double> efficientMeanFilter(CImg<double> input, int filterSize);
CImg<double> efficientGaussianFilter(CImg<double> input, int filterSize);
 
int main(int argc, char **argv) { 
	if (argc < 4) {
		cout << "Insufficent number of arguments." << endl;
		cout << "filters problemID inputfile outputfile" << endl;
		return -1;
	} 
	char* inputFile = argv[2];
	char* outputFile = argv[3];
	cout << "In: " << inputFile << "  Out: " << outputFile << endl; 
	CImg<double> input(inputFile);
	CImg<double> output; 
	if (!strcmp(argv[1], "avgGray")) {
		cout << "# Average Grayscale" << endl;
		if (input.spectrum() != 3) {
			cout << "INPUT ERROR: Input image is not a color image!" << endl;
			return -1;
		}
		output = averageGrayscale(input);
		output.save(outputFile);
	} else if (!strcmp(argv[1], "simpleBW")) {
		cout << "# Simple Threshold Black and White" << endl;
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		output = simpleBW(input);
		output.save(outputFile);
	} else if (!strcmp(argv[1], "advancedBW")) {
		cout << "# Advanced Threshold Black and White" << endl;
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		output = advancedBW(input);
		output.save(outputFile);
	} else if (!strcmp(argv[1], "uniNoise")) {
		cout << "# Uniform Noise" << endl;
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		output = uniformNoise(input);
		output.save(outputFile);
	} else if (!strcmp(argv[1], "gaussNoise")) {
		cout << "# Gaussian Noise" << endl;
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		if (argc != 5) {
			cout << "INPUT ERROR: Provide sigma as additional argument!" << endl;
			return -1;
		}
		output = gaussianNoise(input, atof(argv[4]));
		output.save(outputFile);
	} else if (!strcmp(argv[1], "spNoise")) {
		cout << "# Salt & Pepper Noise" << endl;
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		output = saltAndPepperNoise(input);
		output.save(outputFile);
	} else if (!strcmp(argv[1], "meanFilter")) {
		cout << "# Mean Filter" << endl;
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		if (argc != 5) {
			cout << "INPUT ERROR: Provide filter size as additional argument!" << endl;
			return -1;
		}
		std::clock_t start = std::clock();   	 
		output = meanFilter(input, atoi(argv[4]));
		std::cout << "Time Elapsed: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
		output.save(outputFile);
	} else if (!strcmp(argv[1], "gaussFilter")) {
		cout << "# Gaussian Filter" << endl;
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		if (argc != 5) {
			cout << "INPUT ERROR: Provide sigma as additional argument!" << endl;
			return -1;
		}
		std::clock_t start = std::clock();   	 
		output = gaussianFilter(input, atof(argv[4]));
		std::cout << "Time Elapsed: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
		output.save(outputFile);
	} else if (!strcmp(argv[1], "medianFilter")) {
		cout << "# Median Filter" << endl;
		if (input.spectrum() != 1) {
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		if (argc != 5) {
			cout << "INPUT ERROR: Provide filter size as additional argument!" << endl;
			return -1;
		}
		output = medianFilter(input, atoi(argv[4]));
		output.save(outputFile);
	} else if (!strcmp(argv[1], "analysis")) {
		cout << "# Noise Removal Analysis" << endl;
		CImg<double> original(argv[2]);
		CImg<double> filtered(argv[3]); 
		std::clock_t start = std::clock();
		double MSE = meanSquaredError(original, filtered);
		std::cout << "Time Elapsed: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
		cout << "The Mean Squared Error is: " << MSE << endl;
	} else if (!strcmp(argv[1], "advancedFilter")) {
		cout << "# Separable Kernel Convolutions" << endl; 
		int size = atoi(argv[4]);
		CImg<double> hX(size,1,1,1);
		CImg<double> hY(1,size,1,1);
		cout << "Please give the values of the first 1 - D kernel:"<< endl;
		for(int i = 0; i < size; i++) {
			cout << i + 1 << ": ";
			cin >> hX(i,0,0,0);
		}
		cout << "Please give the values of the second 1 - D kernel:"<< endl;
		for(int i = 0; i < size; i++) {
			cout << i + 1 << ": ";
			cin >> hY(0,i,0,0);
		}	 
		std::clock_t start = std::clock();
		output = separableKernelFilter(input, hX, hY); 
		std::cout << "Time Elapsed: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
		output.save(outputFile);  
	} else if (!strcmp(argv[1], "dynamicFilter")) {
		cout << "# Dynamic Box Filter" << endl; 
		std::clock_t start = std::clock();
		output = efficientMeanFilter(input, atoi(argv[4]));    	 
		std::cout << "Time Elapsed: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
		output.save(outputFile); 
	} else if (!strcmp(argv[1], "fastGaussian")) {
		cout << "# Fast Gaussian Smoothing" << endl; 
		std::clock_t start = std::clock();
		output = efficientGaussianFilter(input, atoi(argv[4]));
		std::cout << "Time Elapsed: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
		output.save(outputFile);
	} else {
		cout << "Unknown input command" << endl;
	} 
	return 0; 
}
 
//Basic image operations
CImg<double> averageGrayscale(CImg<double> input) { 
	//Creates a new grayscale image with same size as the input initialized to all 0s (black)
	CImg<double> output(input.width(), input.height(), 1, 1);
	for (int x = 0; x < input.width(); x++) {
		for (int y = 0; y < input.height(); y++) {
			double R = input(x, y, 0, 0);
			double G = input(x, y, 0, 1);
			double B = input(x, y, 0, 2);
			//Weighted values for R,G and B
			output(x, y, 0, 0) = (0.3 * R + 0.6 * G + 0.1 * B);
		}
	}
	return output; 
}
 
CImg<double> simpleBW(CImg<double> input) { 
	//Creates a new Black and White image with same size as the input initialized to all 0s (black)
	CImg<double> output(input, "xyzc", 0);
	int width = input.width();
	int height = input.height();
	//Compare each pixel to the threshold 127; if larger than the threshold, set to 255, and otherwise set to 0.
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			output(x, y, 0, 0) = input(x, y, 0, 0) < 127 ? 0 : 255;
		}
	}
	return output; 
}
 
CImg<double> advancedBW(CImg<double> input) {
	//Creates a new grayscale image with same size as the input initialized to all 0s (black)
	int error;
	CImg<double> output(input.width(), input.height(), 1, 1);
	int width = input.width();
	int height = input.height();
	//Compare each pixel to the threshold 127; if larger than the threshold, set to 255, and otherwise set to 0.
	//add 0.4e to the pixel to the right of p, 0.3e to the pixel below p, 0.2e to the pixel below and to the left of p,
	//and 0.1e to the pixel below and to the right of p
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			output(x, y, 0, 0) = input(x, y, 0, 0) < 127 ? 0 : 255;
			error=input(x,y,0,0)-output(x,y,0,0);
			//changes to the right pixel
			if(x+1<width)                        	 
				input(x+1,y,0)=input(x+1,y,0)+0.4*error;
			//changes to the bottom pixel
			if(y+1<height)                       	 
				input(x,y+1,0)=input(x,y+1,0)+0.3*error;
			//changes to the pixel to the right and below
			if(x+1<width && y+1<height)  	 
				input(x+1,y+1,0)=input(x+1,y+1,0)+0.1*error;
			//changes to the pixel to the left and below
			if(x>0 && y+1<height)                 	 
				input(x-1,y+1,0)=input(x-1,y+1,0)+0.2*error;
		}
	}
	return output;
}
 
//Adding noise
CImg<double> uniformNoise(CImg<double> input) {
	double noise;
	int randomNum,sum;
	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input, "xyzc", 0);
	//Adding noise uniformly sampled at random number between -20 and +20(given)
	for (int x=0;x<input.width();x++) {
		for(int y=0;y<input.height();y++) {
			//Generate random number between -20 to 20  		 
			randomNum= (rand()%40)-20;
			output(x,y,0,0) = input(x,y,0,0)+randomNum;
			//clamping the noisy values at 0 and 255, so that they don't become negative or larger than 255
			if(output(x,y,0,0)>255)
				output(x,y,0,0) =255;
			if(output(x,y,0,0)<0)
				output(x,y,0,0)=0;
		}
	}
	return output;
}
 
CImg<double> gaussianNoise(CImg<double> input, double sigma) {
	double PI = 3.14159265359,noiseGaussian;
	double v1, v2, s, x;
	int stage = 0;
	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input, "xyzc", 0); 
	//Adding noise s from a zero-mean Gaussian with a given standard deviation
	for(int i=0;i<output.width();i++){
		for(int j=0;j<input.height();j++){
			if(stage == 0) {
				do {
					double u1 = (double)rand() / RAND_MAX;
					double u2 = (double)rand() / RAND_MAX;
					v1 = 2 * u1 - 1;
					v2 = 2 * u2 - 1;
					s = v1 * v1 + v2 * v2;
				} while(s >= 1 || s == 0);			
				x = v1 * sqrt(-2 * log(s) / s);
			} 
			else
				x = v2 * sqrt(-2 * log(s) / s);
			//Modify stage by subtracting by 1
			stage = 1 - stage;
			//The additive value for Gaussian noise
			noiseGaussian=x*sigma;
			output(i,j,0,0)=input(i,j,0,0)+noiseGaussian;
			//clamping the noisy values at 0 and 255, so that they don't become negative or larger than 255
			if(output(i,j,0,0)>255)
				output(i,j,0,0)=255;
			if(output(i,j,0,0)<0)
				output(i,j,0,0)=0;
		}
	}
	return output;
}
 
 
CImg<double> saltAndPepperNoise(CImg<double> input) { 
	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input, "xyzc", 0);
	/*For a given pixel p, this type of noise rst decides whether to add noise (with some small probability P) or whether to keep the pixel value
	*as-is. If it decides to add noise, then it sets the pixel value to either 0 (with probability 50%) or 255 (with probability 50%). This mimicks
	*pixel values that are \stuck" at 0 or 255 due to, for example, defective pixels in the camera sensor.
	*/
	int rp,bw;
	for (int x=0;x<input.width();x++){
		for(int y=0;y<input.height();y++){
			output(x,y,0,0) = input(x,y,0,0);
			//Finding probability(Low)
			//Assumption: : Random number between 1-10 to achieve a probability of (0.1)
			rp = (rand()%10)+1;
			//check if random pixel is 5
			if(rp == 5) {
				bw = (rand()%2)+1;  
				if(bw == 1)    	
					output(x,y,0,0) = 0;
				else
					output(x,y,0,0) = 255;
			}
		}
	}
	return output;
}
 
//Filtering
CImg<double> filter(CImg<double> input, CImg<double> filter){
	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input.width(), input.height(), 1, 1);
	//Convolve with filter
	int width = input.width();
	int height = input.height();
	int filterSize = filter.width();
	int k = (filterSize - 1) / 2;	
	for(int i = k; i < width - k; i++){
		for(int j = k; j < height - k; j++){
			output(i,j,0,0) = 0;
			for(int x = 0; x < filterSize; x++) {
				for(int y = 0; y < filterSize; y++) {
					output(i,j,0,0) += input(i - k + x,j - k + y,0,0) * filter(x,y,0,0);
				} 
			}
		}
	}
	return output;
}

CImg<double> meanFilter(CImg<double> input, int filterSize){
	//Creates a new grayscale image (just a matrix) to be our filter
	CImg<double> H(filterSize, filterSize, 1, 1); 
	//Fill filter values
	for(int i = 0; i < filterSize; i++) {
		for(int j = 0; j < filterSize; j++) {
			H(i,j,0,0) = (double) 1 / (filterSize * filterSize);
		}
	}	
	//Convole with filter and return
	return filter(input, H);	
}

CImg<double> gaussianFilter(CImg<double> input, double sigma){
	//Determine filter size, see Part 4 2b for how
	int filterSize = (int) round(3 * sigma);
	int k = (filterSize - 1) / 2;		
	double pi = 3.14159265359;	
	//Creates a new grayscale image (just a matrix) to be our filter 
	CImg<double> H(filterSize, filterSize, 1, 1);	
	double total = 0;
	//Fill filter values	
	for(int i = 0; i < filterSize; i++) {
		double x = i - k;
		for(int j = 0; j < filterSize; j++) {
			double y = j - k;
			double power = (double) -(x * x + y * y) / (2 * sigma * sigma);
			H(i,j,0,0) = (double) exp(power) / (2 * pi * sigma * sigma);
		}
	}	
	//Convole with filter and return
	return filter(input, H);
}

CImg<double> medianFilter(CImg<double> input, int size){
	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input, "xyzc", 0);
	int width = input.width();
	int height = input.height();
	int k = (size - 1) / 2;
	double a[size * size];
	//Do median filtering
	for (int i = 0; i < input.width(); i++) {
		for (int j = 0; j < input.height(); j++) {
			for (int x = 0; x < size; x++) {
				for (int y = 0; y < size; y++) {
					if (i - k + x >= 0 && j - k + y >= 0 && i - k + x < width && j - k + y < height) {
						a[x * size + y] = input(i - k + x, j - k + y, 0, 0);
					}
					else {
						a[x * size + y] = 0;
					}
				}
			}
			sort(a, a + size * size);
			output(i, j, 0, 0) = a[(size*size - 1) / 2];
		}
	}
	return output;
}

double meanSquaredError(CImg<double> original, CImg<double> filtered) {
	double MSE = 0;
	int width = original.width();
	int height = original.height();
	for(int i = 0; i < width; i++){
		for(int j = 0; j < height; j++){			 
			double diff = original(i,j,0,0) - filtered(i,j,0,0);
			MSE += (diff * diff) / (width * height);
		}
	}
	return MSE;
}

CImg<double> separableKernelFilter(CImg<double> input, CImg<double>  filterX, CImg<double>  filterY){
	int filterSize = filterX.width();
	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> firstPass(input.width(), input.height(), 1, 1);
	CImg<double> secondPass(input.width(), input.height(), 1, 1);
	//Convolve with filter
	int width = input.width();
	int height = input.height();
	int k = (filterSize - 1) / 2;	
	for(int i = k; i < width - k; i++){
		for(int j = k; j < height - k; j++){
			firstPass(i,j,0,0) = 0;
			for(int x = 0; x < filterSize; x++) {
				firstPass(i,j,0,0) += input(i,j - k + x,0,0) * filterY(0,x,0,0);
			}
		}
	}	
	for(int j = k; j < height - k; j++){
		for(int i = k; i < width - k; i++){
			secondPass(i,j,0,0) = 0;
			for(int x = 0; x < filterSize; x++) {
				secondPass(i,j,0,0) += firstPass(i - k + x,j,0,0) * filterX(x,0,0,0);
			}
		}
	}
	return secondPass;
}

CImg<double> efficientMeanFilter(CImg<double> input, int filterSize){	
	double total = (double) 1 / filterSize;	
	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> firstPass(input.width(), input.height(), 1, 1);
	CImg<double> secondPass(input.width(), input.height(), 1, 1);
	//Convolve with filter
	int width = input.width();
	int height = input.height();
	int k = (filterSize - 1) / 2;	
	for(int i = k; i < width - k; i++){
		for(int j = k; j < height - k; j++){			 
			if(j > k) {
				firstPass(i,j,0,0) = firstPass(i,j - 1,0,0) - (input(i,j - 1 - k,0,0) * total) + (input(i,j + k,0,0) * total);
			} else {
				firstPass(i,j,0,0) = 0;				
				for(int x = 0; x < filterSize; x++) {
					firstPass(i,j,0,0) += input(i,j - k + x,0,0) * total;
				}
			}
		}
	}	
	for(int j = k; j < height - k; j++){
		for(int i = k; i < width - k; i++){
			if(i > k) {
				secondPass(i,j,0,0) = secondPass(i - 1,j,0,0) - (firstPass(i - 1 - k,j,0,0) * total) + (firstPass(i + k,j,0,0) * total);								
			} else {
				secondPass(i,j,0,0) = 0;
				for(int x = 0; x < filterSize; x++) {
					secondPass(i,j,0,0) += firstPass(i - k + x,j,0,0) * total;
				}
			}
		}
	}
	return secondPass;	
}

CImg<double> efficientGaussianFilter(CImg<double> input, int filterSize){	
	// The filter is applied for 3 iterations by default
	int repetitions = 3;
	std::cout << "Number of iterations the filter should be applied (3 by default): ";
	std::string in;
	std::getline(std::cin, in);
	if (!in.empty()) {
		std::istringstream stream(in);
		stream >> repetitions;
	}	
	CImg<double> output = efficientMeanFilter(input, filterSize);	
	for(int i = 0; i < repetitions - 1; i++) {
		output = efficientMeanFilter(output, filterSize);
	}	
	return output;	
}
