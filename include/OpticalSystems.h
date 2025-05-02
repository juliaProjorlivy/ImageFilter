#pragma once
#include "FourierTransform.h"
#include "LightSources.h"
#include <fftw3.h>
#include <cmath>
#include <cstring>

class OpticalSystemBase
{
protected:
    int width_;
    int height_;
    double R_;
    double Lambda_;

    FourieTransform fourie_;
public:
    OpticalSystemBase(int width, int height, double R, double Lambda) : width_(width), height_(height), R_(R), Lambda_(Lambda), fourie_(width, height) {}
    OpticalSystemBase(const FourieTransform &fourie, double R, double Lambda) : width_(fourie.width_), height_(fourie.height_), R_(R), Lambda_(Lambda), fourie_(fourie.width_, fourie.height_) 
    {
        std::memcpy(fourie_.in_, fourie.in_, sizeof(fftw_complex)*width_*height_);
        std::memcpy(fourie_.out_, fourie.out_, sizeof(fftw_complex)*width_*height_);
        // // Copy data manually
        // for (int i = 0; i < height_; ++i) {
        //     for(int j = 0; j < width_; ++j) {
        //         fourie_.in_[i * width_ + j][0]  = fourie.in_[i * width_ + j][0];
        //         fourie_.in_[i * width_ + j][1]  = fourie.in_[i * width_ + j][1];
        //         fourie_.out_[i * width_ + j][0] = fourie.out_[i *width_ + j][0];
        //         fourie_.out_[i * width_ + j][1] = fourie.out_[i *width_ + j][1];
        //     }
        // }
    }

    virtual ~OpticalSystemBase() = default;

    void changeR(double R) 
    {
        R_ = std::clamp(R + R_, 1.0, static_cast<double>(std::min(width_, height_)));
    }

    void changeLambda(double lambda) 
    {
        Lambda_ += lambda;
    }

    // Функция для описания дифракции
    virtual void apply() = 0;

    // Вычисление прямого фурье преобразования квадрата функции рассеяния точки
    void computePSFSquaredFT() 
    {
        for (int i = 0; i < width_ * height_; ++i) 
        {
            fourie_.in_[i][0] = fourie_.out_[i][0] * fourie_.out_[i][0] - fourie_.out_[i][1] * fourie_.out_[i][1];
            fourie_.in_[i][1] = fourie_.out_[i][0] * fourie_.out_[i][1] + fourie_.out_[i][1] * fourie_.out_[i][0];
        }

        fourie_.transform(FFTW_FORWARD);
    }

    const fftw_complex *computeIntencity(const LightSource &source) 
    {
        apply();
        computePSFSquaredFT();

        for (int i = 0; i < width_ * height_; ++i) 
        {
            double source_re = source.getOut(i, 0);
            double source_im = source.getOut(i, 1);
            fourie_.in_[i][0] = fourie_.out_[i][0] * source_re - fourie_.out_[i][1] * source_im;
            fourie_.in_[i][1] = fourie_.out_[i][0] * source_im + fourie_.out_[i][1] * source_re;
        }

        fourie_.transform(FFTW_BACKWARD);
        return fourie_.out_;
    }

    double getR() const 
    {
        return R_; 
    }

    double getLambda() const 
    {
        return Lambda_;
    }

    int getSize() const 
    {
        return width_; 
    }
};

class CircleAperture : public OpticalSystemBase 
{
public:
    CircleAperture(int width, int height, double R, double Lambda) : OpticalSystemBase(width, height, R, Lambda) {}
    CircleAperture(const FourieTransform &fourie, double R, double Lambda) : OpticalSystemBase(fourie, R, Lambda) {}
    ~CircleAperture() override = default;

    // Применение круглой диафрагмы (обнуление пикселей за пределами радиуса R_)
    void apply() override 
    {
        for (int i = 0; i < height_; ++i) 
        {
            double y = i - height_ / 2.0;

            for (int j = 0; j < width_; ++j) 
            {
                double x = j - width_ / 2.0;
                double r = sqrt(x * x + y * y);

                // Super-Gaussian aperture (sharp but smooth)
                fourie_.in_[i * width_ + j][0] = exp(-pow(r / R_, 8.0));
                // Hann window (reduces edge artifacts)
                double window = 0.5 * (1 + cos(M_PI * r / (std::max(width_, height_) / 2.)));
                fourie_.in_[i * width_ + j][0] *= window;
                fourie_.in_[i * width_ + j][1] = 0.0;
            }
        }

        fourie_.transform(FFTW_BACKWARD);
    }

};

class SerfaceEncoding : public OpticalSystemBase 
{
public:
    SerfaceEncoding(int width, int height, double R, double Lambda) : OpticalSystemBase(width, height, R, Lambda) {}
    ~SerfaceEncoding() override = default;

    // Применение фазового кодирования (g(x,y) = sign(x)*|x|^2.5 + sign(y)*|y|^2.5)
    void apply() override 
    {
        double k = 2 * M_PI / Lambda_;
        for (int i = 0; i < height_; ++i) 
        {
            double y = i - height_ / 2.0;

            for (int j = 0; j < width_; ++j) 
            {
                double x = j - width_ / 2.0;
                double r = sqrt(x * x + y * y);

                // Super-Gaussian aperture (sharp but smooth)
                double val = exp(-pow(r / R_, 8.0));
                // Hann window (reduces edge artifacts)
                double window = 0.5 * (1 + cos(M_PI * r / (std::max(width_, height_) / 2.)));
                // Фазовая модуляция
                double phase = copysign(pow(std::fabs(x), 2.5), x) + 
                    copysign(pow(std::fabs(y), 2.5), y);

                fourie_.in_[i * width_ + j][0] = val * window * cos(k * phase);
                fourie_.in_[i * width_ + j][1] = val * window * sin(k * phase);
            }
        }

        fourie_.transform(FFTW_BACKWARD);
    }
};