#pragma once
#include "OpticalSystems.h"
#include <cmath>
#include <fftw3.h>
#include <algorithm>

class InverseFiltering
{
private:
    int width_;
    int height_;
    double Lambda_;
    FourieTransform fourie_;

    const double epsilon_ = 1e-6;
public:
    InverseFiltering(int width, int height, double Lambda) : width_(width), height_(height), Lambda_(Lambda), fourie_(width, height) {}
    ~InverseFiltering() = default;

    FourieTransform &&filter(const fftw_complex *distorted_image, const OpticalSystemBase &optical_system, const LightSource &source) {

        // 1. Получение ~B * ~h
        for (int i = 0; i < width_ * height_; ++i)
        {
            fourie_.in_[i][0] = distorted_image[i][0]; // Re
            fourie_.in_[i][1] = distorted_image[i][0]; // Im
        }

        fourie_.transform(FFTW_FORWARD);

        // 2. Получение функции рассеяния точки
        for (int i = 0; i < width_ * height_; ++i)
        {
            double br_re = source.getOut(i, 0);
            double br_im = source.getOut(i, 1);
            double denom = br_re * br_re + br_im * br_im + epsilon_;

            fourie_.in_[i][0] = (fourie_.out_[i][0] * br_re + fourie_.out_[i][1] * br_im) / denom; // Re
            fourie_.in_[i][1] = (fourie_.out_[i][1] * br_re - fourie_.out_[i][0] * br_im) / denom; // Im
        }

        fourie_.transform(FFTW_BACKWARD);

        // 3. Получение зрачковой функции (теперь ее нужно отфильтровать от поверхности g(x,y) = sign(x)*|x|^2.5 + sign(y)*|y|^2.5)
        for (int i = 0; i < width_ * height_; ++i)
        {
            double re = fourie_.out_[i][0];
            double im = fourie_.out_[i][1];
            double denom = re * re + im * im + epsilon_;

            fourie_.in_[i][0] = (re * re + im * im) / denom; // Re
            fourie_.in_[i][1] = (im * re - re * im) / denom; // Im
        }

        fourie_.transform(FFTW_FORWARD);

        double k = 2 * M_PI / Lambda_;
        
        for (int i = 0; i < height_; ++i)
        {
            double y = i - height_ / 2.0;
            
            for (int j = 0; j < width_; ++j)
            {
                double x = j - width_ / 2.0;

                // Фазовая модуляция
                double phase = copysign(pow(std::fabs(x), 2.5), x) + copysign(pow(std::fabs(y), 2.5), y);

                fourie_.in_[i * width_ + j][0] = fourie_.out_[i * width_ + j][0] / (cos(k * phase) + epsilon_);
                fourie_.in_[i * width_ + j][1] = fourie_.out_[i * width_ + j][1] / (sin(k * phase) + epsilon_) - fourie_.in_[i * width_ + j][0];
            }
        }

        fourie_.transform(FFTW_BACKWARD);

        return std::move(fourie_);
    }
};