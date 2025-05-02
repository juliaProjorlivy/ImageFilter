#pragma once
#include "FourierTransform.h"
#include <vector>
#include <fftw3.h>
#include <algorithm>
#include <cstring>

class LightSource 
{
protected:
    int width_;
    int height_;
    int lenX_;
    int lenY_;
    std::vector<double> brightness_;
    std::vector<double> brightness_squared_;
    FourieTransform fourie_;

public:
    LightSource(int width, int height) : width_(width), height_(height), brightness_(width * height, 0.0), brightness_squared_(width * height, 0.0), fourie_(width, height) {}

    virtual ~LightSource() = default;

    virtual void setBrightness() = 0;

    void setSquaredBrightness() 
    {
        for (int i = 0; i < width_ * height_; ++i) 
        {
            brightness_squared_[i] = brightness_[i] * brightness_[i];
        }
    }

    void setFourieBrightnessSquared() 
    {
        // Копируем B² в комплексный массив для FFT
        for (int i = 0; i < width_ * height_; ++i) 
        {
            fourie_.in_[i][0] = brightness_squared_[i]; // Re
            fourie_.in_[i][1] = 0.0;                    // Im
        }

        // Вычисляем Фурье-образ B²
        fourie_.transform(FFTW_FORWARD);
    }

    void changeLenX(int x)
    {
        lenX_ = std::clamp(lenX_ + x, 1, width_);
        setBrightness();
        setSquaredBrightness();
        setFourieBrightnessSquared();
    }

    void changeLenY(int y)
    {
        lenY_ = std::clamp(lenY_ + y, 1, width_);
        setBrightness();
        setSquaredBrightness();
        setFourieBrightnessSquared();
    }

    int getLenX()
    {
        return lenX_;
    }
    
    int getLenY()
    {
        return lenY_;
    }

    double getOut(int index, int sign) const 
    {
        return fourie_.out_[index][sign];
    }

    const fftw_complex* getFourieData() const 
    {
        return fourie_.out_;
    }
};

class PointSource : public LightSource 
{
public:

    PointSource(int width, int height) : LightSource(width, height) 
    {
        setBrightness();
        setSquaredBrightness();
        // F[delta] = 1
        setFourieBrightnessSquared();
    }
    ~PointSource() = default;

    void setBrightness() override 
    {
        std::fill(brightness_.begin(), brightness_.end(), 0.0);
        brightness_[(height_ / 2) * width_ + (width_ / 2)] = 1.0; // Центр изображения
    }
};

class LineSource : public LightSource 
{
public:
    LineSource(int width, int height, int lenX) : LightSource(width, height)
    {
        lenX_ = lenX;
        setBrightness();
        setSquaredBrightness();
        setFourieBrightnessSquared();
    }

    void setlenX(int lenX) 
    {
        lenX_ = std::clamp(lenX, 1, width_);
    }

    void setBrightness() override
    {
        std::fill(brightness_.begin(), brightness_.end(), 0.0);

        int x_left = (width_ - lenX_) / 2;
        int x_right = width_ - x_left;

        if (width_ < lenX_)
        {
            x_left = 0; 
            x_right = width_;
        }

        for (int x = x_left; x < x_right; ++x)
        {
            brightness_[(height_ / 2) * width_ + x] = 1.0; // Средняя строка
        }
    }
};

class CrossSource : public LightSource
{
public:
    CrossSource(int width, int height, int lenX, int lenY) : LightSource(width, height)
    {
        lenX_ = lenX;
        lenY_ = lenY;
        setBrightness();
        setSquaredBrightness();
        setFourieBrightnessSquared();
    }

    void setlenX(int lenX) 
    {
        lenX_ = std::clamp(lenX, 1, width_);
    }
    
    void setlenY(int lenY) 
    {
        lenY_ = std::clamp(lenY, 1, height_);
    }

    void setBrightness() override
    {
        std::fill(brightness_.begin(), brightness_.end(), 0.0);

        // Горизонтальная линия
        int x_left = (width_ - lenX_) / 2;
        int x_right = width_ - x_left;
        
        for (int x = x_left; x < x_right; ++x)
        {
            brightness_[(height_ / 2) * width_ + x] = 1.0;
        }

        // Вертикальная линия
        int y_high = (height_ - lenY_) / 2;
        int y_low = height_ - y_high;
        
        for (int y = y_low; y < y_high; ++y)
        {
            brightness_[y * width_ + (width_ / 2)] = 1.0;
        }
    }
};