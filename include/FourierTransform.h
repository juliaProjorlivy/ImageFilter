#pragma once
#include <fftw3.h>
#include <algorithm>

struct FourieTransform 
{
    int width_;
    int height_;
    fftw_complex* in_;
    fftw_complex* out_;

    FourieTransform(int width, int height): width_(width), height_(height) 
    {
        in_ = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * width * height);
        out_ = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * width * height);
    }
    
    ~FourieTransform()
    { 
        fftw_free(in_); 
        fftw_free(out_);
    }

    // Сдвиг квадрантов для корректного отображения FFT
    void fftshift(fftw_complex* data) 
    {
        int half_w = width_ / 2;
        int half_h = height_ / 2;

        for (int y = 0; y < half_h; ++y) 
        {
            for (int x = 0; x < half_w; ++x) 
            {
                // Индексы для четырех квадрантов
                int src_idx1 = y * width_ + x;
                int trg_idx1 = (y + half_h) * width_ + (x + half_w);

                int src_idx2 = y * width_ + (x + half_w);
                int trg_idx2 = (y + half_h) * width_ + x;

                // Обмен данными
                std::swap(data[src_idx1][0], data[trg_idx1][0]);
                std::swap(data[src_idx1][1], data[trg_idx1][1]);

                std::swap(data[src_idx2][0], data[trg_idx2][0]);
                std::swap(data[src_idx2][1], data[trg_idx2][1]);
            }
        }
    }

    void transform(int sign = FFTW_FORWARD) 
    {
        fftw_plan plan = fftw_plan_dft_2d(width_, height_, in_, out_, sign, FFTW_ESTIMATE);
        fftw_execute(plan);   // Результат в 'out'
        fftshift(out_);    // Центрируем частоты
    }
};