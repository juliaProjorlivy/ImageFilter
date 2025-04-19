#ifndef IMAGE_H
#define IMAGE_H
#include <SFML/Window/Keyboard.hpp>
#include <fftw3.h>
#include <SFML/Graphics.hpp>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string.h>

class DotObject
{
private:
    int nx_; // Number of pixels_
    int ny_;
    double R_;
    fftw_complex *in_;
    fftw_complex *out_;

    sf::Uint8* pixels_;
    fftw_plan plan_;

public:
    DotObject(int nx, int ny, double R) : nx_(nx), ny_(ny), R_(R)
    {
        in_  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx_ * ny_);
        out_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx_ * ny_);
        pixels_ = new sf::Uint8[nx_ * ny_ * 4]; // RGBA array
    }

    DotObject(DotObject& dot) : nx_(dot.nx_), ny_(dot.ny_), R_(dot.R_)
    {
        in_  = (fftw_complex*) memcpy(in_, dot.in_, sizeof(fftw_complex) * nx_ * ny_);
        out_ = (fftw_complex*) memcpy(out_, dot.out_, sizeof(fftw_complex) * nx_ * ny_);
        pixels_ = (sf::Uint8*) memcpy(pixels_, dot.pixels_, sizeof(sf::Uint8) * nx_ * ny_ * 4);
        plan_ = dot.plan_;
    }

    DotObject(DotObject&& dot) : nx_(dot.nx_), ny_(dot.ny_), R_(dot.R_), in_(dot.in_), out_(dot.out_), pixels_(dot.pixels_), plan_(dot.plan_)
    {
        dot.in_ = nullptr;
        dot.out_ = nullptr;
        dot.pixels_ = nullptr;
        dot.plan_ = nullptr;
    }

    DotObject& operator=(const DotObject& dot) = delete;
    DotObject& operator=(const DotObject&& dot) = delete;

    ~DotObject()
    {
        fftw_destroy_plan(plan_);
        fftw_free(in_);
        fftw_free(out_);
        delete[] pixels_;
    }


    void updateImage(double R, sf::Texture& texture)
    {
        // Updating the radius
        R_ = R;
        // Calculating the new in_ and out_
        transmission_func();

        fftshift(in_, nx_, ny_); // Shift before FFT
        plan_ = fftw_plan_dft_2d(nx_, ny_, in_, out_, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_);
        fftshift(out_, nx_, ny_); // Shift before FFT

        fill_image();

        // TODO: probably should remove graphics from class (pixels_ aslo should be removed)
        texture.setSmooth(true); // Enable smooth display interpolation
        texture.update(pixels_);
    }

private:
    // The specific Fourier Transform of FFTW
    void fftshift(fftw_complex* data, int nx, int ny)
    {
        // Swap quadrants for rectangular arrays
        int half_nx = nx / 2;
        int half_ny = ny / 2;
        for (int i = 0; i < half_nx; ++i)
        {
            for (int j = 0; j < half_ny; ++j)
            {
                // Swap quadrant 1 and 3
                std::swap(data[i * ny + j], data[(i + half_nx) * ny + (j + half_ny)]);
                // Swap quadrant 2 and 4
                std::swap(data[i * ny + (j + half_ny)], data[(i + half_nx) * ny + j]);
            }
        }
    }

    virtual void transmission_func()
    {
        // Soft-edged aperture + Hann window
        for (int i = 0; i < ny_; ++i)
        {
            for (int j = 0; j < nx_; ++j)
            {
                double x = i - nx_ / 2.0;
                double y = j - ny_ / 2.0;
                double r = sqrt(x * x + y * y);

                // Super-Gaussian aperture (sharp but smooth)
                in_[i * ny_ + j][0] = exp(-pow(r / R_, 8.0));
                // Hann window (reduces edge artifacts)
                double window = 0.5 * (1 + cos(M_PI * r / (std::max(nx_, ny_) / 2.)));
                in_[i * ny_ + j][0] *= window;
                in_[i * ny_ + j][1] = 0.0;
            }
        }
    }

    virtual void fill_image()
    {
        double max_log_intensity = 0.0;
        for (int i = 0; i < ny_; ++i)
        {
            for (int j = 0; j < nx_; ++j)
            {
                double real = out_[i * ny_ + j][0];
                double imag = out_[i * ny_ + j][1];
                double log_intensity = log(1 + 100 * (real * real + imag * imag));
                if (log_intensity > max_log_intensity) max_log_intensity = log_intensity;
            }
        }

        // Store result with smooth interpolation
        for (int i = 0; i < ny_; ++i)
        {
            for (int j = 0; j < nx_; ++j)
            {
                double real = out_[i * ny_ + j][0];
                double imag = out_[i * ny_ + j][1];
                double log_intensity = log(1 + 100 * (real * real + imag * imag));
                sf::Uint8 value = static_cast<sf::Uint8>(255 * log_intensity / max_log_intensity);
                pixels_[(i * ny_ + j) * 4]     = value;
                pixels_[(i * ny_ + j) * 4 + 1] = value;
                pixels_[(i * ny_ + j) * 4 + 2] = value;
                pixels_[(i * ny_ + j) * 4 + 3] = 255;
            }
        }
    }
};


#endif // !IMAGE_H
