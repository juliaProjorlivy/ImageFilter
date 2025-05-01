#include <SFML/Graphics.hpp>
#include <fftw3.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>

// ==================================================
// Класс Фурье преобразования
// ==================================================
struct FourieTransform {
    int width_;
    int height_;
    fftw_complex* in_;
    fftw_complex* out_;

    FourieTransform(int width, int height): width_(width), height_(height) {
        in_ = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * width * height);
        out_ = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * width * height);
    }
    ~FourieTransform(){ fftw_free(in_); fftw_free(out_);}

    // Сдвиг квадрантов для корректного отображения FFT
    void fftshift(fftw_complex* data) {
        int half_w = width_ / 2;
        int half_h = height_ / 2;

        for (int y = 0; y < half_h; ++y) {
            for (int x = 0; x < half_w; ++x) {
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

    void transform(int sign = FFTW_FORWARD) {
        fftw_plan plan = fftw_plan_dft_2d(width_, height_, in_, out_, sign, FFTW_ESTIMATE);
        fftw_execute(plan);   // Результат в 'out'
        fftshift(out_);    // Центрируем частоты
    }
};

// ==================================================
// Базовый класс для источников света
// ==================================================
class LightSource {
protected:
    int width_;
    int height_;
    std::vector<double> brightness_;
    std::vector<double> brightness_squared_;
    FourieTransform fourie_;

public:
    LightSource(int width, int height) : width_(width), height_(height), brightness_(width * height, 0.0), brightness_squared_(width * height, 0.0), fourie_(width, height) {}

    virtual ~LightSource() = default;

    virtual void setBrightness() = 0;

    void setSquaredBrightness() {
        for (int i = 0; i < width_ * height_; ++i) {
            brightness_squared_[i] = brightness_[i] * brightness_[i];
        }
    }

    void setFourieBrightnessSquared() {
        // Копируем B² в комплексный массив для FFT
        for (int i = 0; i < width_ * height_; ++i) {
            fourie_.in_[i][0] = brightness_squared_[i]; // Re
            fourie_.in_[i][1] = 0.0;                    // Im
        }
        // Вычисляем Фурье-образ B²
        fourie_.transform(FFTW_FORWARD);
    }

    double getOut(int index, int sign) const {return fourie_.out_[index][sign];}
};

// ==================================================
// Точечный источник (дельта-функция в центре)
// ==================================================
class PointSource : public LightSource {
public:

    PointSource(int width, int height) : LightSource(width, height) {
        setBrightness();
        setSquaredBrightness();
        // F[delta] = 1
        setFourieBrightnessSquared();
    }
    ~PointSource() = default;

    void setBrightness() override {
        std::fill(brightness_.begin(), brightness_.end(), 0.0);
        brightness_[(height_ / 2) * width_ + (width_ / 2)] = 1.0; // Центр изображения
    }
};

// ==================================================
// Линия вдоль горизонтальной оси (X)
// ==================================================
class LineSource : public LightSource {
    int lenX_;
public:
    LineSource(int width, int height, int lenX) : LightSource(width, height), lenX_(lenX) {
        setBrightness();
        setSquaredBrightness();
        setFourieBrightnessSquared();
    }

    void setlenX(int lenX) {lenX_ = std::clamp(lenX, 1, width_);}

    void setBrightness() override {
        std::fill(brightness_.begin(), brightness_.end(), 0.0);

        int x_left = (width_ - lenX_) / 2;
        int x_right = width_ - x_left;
        if(width_ < lenX_) {x_left = 0; x_right = width_;}
        for (int x = x_left; x < x_right; ++x) {
            brightness_[(height_ / 2) * width_ + x] = 1.0; // Средняя строка
        }
    }
};

// ==================================================
// Крест (горизонтальная + вертикальная линии)
// ==================================================
class CrossSource : public LightSource {
    int lenX_;
    int lenY_;
public:
    CrossSource(int width, int height, int lenX, int lenY) : LightSource(width, height), lenX_(lenX), lenY_(lenY) {
        setBrightness();
        setSquaredBrightness();
        setFourieBrightnessSquared();
    }

    void setlenX(int lenX) {lenX_ = std::clamp(lenX, 1, width_);}
    void setlenY(int lenY) {lenY_ = std::clamp(lenY, 1, height_);}

    void setBrightness() override {
        std::fill(brightness_.begin(), brightness_.end(), 0.0);

        // Горизонтальная линия
        int x_left = (width_ - lenX_) / 2;
        int x_right = width_ - x_left;
        for (int x = x_left; x < x_right; ++x) {
            brightness_[(height_ / 2) * width_ + x] = 1.0;
        }
        // Вертикальная линия
        int y_high = (height_ - lenY_) / 2;
        int y_low = height_ - y_high;
        for (int y = y_low; y < y_high; ++y) {
            brightness_[y * width_ + (width_ / 2)] = 1.0;
        }
    }
};

class OpticalSystemBase{
protected:
    int width_;
    int height_;
    double R_;
    double Lambda_;

    FourieTransform fourie_;
public:
    OpticalSystemBase(int width, int height, double R, double Lambda) : width_(width), height_(height), R_(R), Lambda_(Lambda), fourie_(width, height) {}
    OpticalSystemBase(const FourieTransform &fourie, double R, double Lambda) : width_(fourie.width_), height_(fourie.height_), R_(R), Lambda_(Lambda), fourie_(fourie.width_, fourie.height_) {
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

    void changeR(double R) {R_ = std::clamp(R + R_, 1.0, static_cast<double>(std::min(width_, height_)));}
    void changeLambda(double lambda) {Lambda_ += lambda;}

    // Функция для описания дифракции
    virtual void apply() = 0;

    // Вычисление прямого фурье преобразования квадрата функции рассеяния точки
    void computePSFSquaredFT() {
        for (int i = 0; i < width_ * height_; ++i) {
            fourie_.in_[i][0] = fourie_.out_[i][0] * fourie_.out_[i][0] - fourie_.out_[i][1] * fourie_.out_[i][1];
            fourie_.in_[i][1] = fourie_.out_[i][0] * fourie_.out_[i][1] + fourie_.out_[i][1] * fourie_.out_[i][0];
        }

        fourie_.transform(FFTW_FORWARD);
    }

    const fftw_complex *computeIntencity(const LightSource &source) {
        apply();
        computePSFSquaredFT();

        for (int i = 0; i < width_ * height_; ++i) {
            double source_re = source.getOut(i, 0);
            double source_im = source.getOut(i, 1);
            fourie_.in_[i][0] = fourie_.out_[i][0] * source_re - fourie_.out_[i][1] * source_im;
            fourie_.in_[i][1] = fourie_.out_[i][0] * source_im + fourie_.out_[i][1] * source_re;
        }

        fourie_.transform(FFTW_BACKWARD);
        return fourie_.out_;
    }
};

class CircleAperture : public OpticalSystemBase {
public:
    CircleAperture(int width, int height, double R, double Lambda) : OpticalSystemBase(width, height, R, Lambda) {}
    CircleAperture(const FourieTransform &fourie, double R, double Lambda) : OpticalSystemBase(fourie, R, Lambda) {}
    ~CircleAperture() override = default;

    // Применение круглой диафрагмы (обнуление пикселей за пределами радиуса R_)
    void apply() override {
        for (int i = 0; i < height_; ++i) {
            double y = i - height_ / 2.0;
            for (int j = 0; j < width_; ++j) {
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

class SerfaceEncoding : public OpticalSystemBase {
public:
    SerfaceEncoding(int width, int height, double R, double Lambda) : OpticalSystemBase(width, height, R, Lambda) {}
    ~SerfaceEncoding() override = default;

    // Применение фазового кодирования (g(x,y) = sign(x)*|x|^2.5 + sign(y)*|y|^2.5)
    void apply() override {
        double k = 2 * M_PI / Lambda_;
        for (int i = 0; i < height_; ++i) {
            double y = i - height_ / 2.0;
            for (int j = 0; j < width_; ++j) {
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

class InverseFiltering {
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
        for (int i = 0; i < width_ * height_; ++i) {
            fourie_.in_[i][0] = distorted_image[i][0]; // Re
            fourie_.in_[i][1] = distorted_image[i][0]; // Im
        }
        fourie_.transform(FFTW_FORWARD);

        // 2. Получение функции рассеяния точки
        for (int i = 0; i < width_ * height_; ++i) {
            double br_re = source.getOut(i, 0);
            double br_im = source.getOut(i, 1);
            double denom = br_re * br_re + br_im * br_im + epsilon_;
            fourie_.in_[i][0] = (fourie_.out_[i][0] * br_re + fourie_.out_[i][1] * br_im) / denom; // Re
            fourie_.in_[i][1] = (fourie_.out_[i][1] * br_re - fourie_.out_[i][0] * br_im) / denom; // Im
        }
        fourie_.transform(FFTW_BACKWARD);

        // 3. Получение зрачковой функции (теперь ее нужно отфильтровать от поверхности g(x,y) = sign(x)*|x|^2.5 + sign(y)*|y|^2.5)
        for (int i = 0; i < width_ * height_; ++i) {
            double re = fourie_.out_[i][0];
            double im = fourie_.out_[i][1];
            double denom = re * re + im * im + epsilon_;
            fourie_.in_[i][0] = (re * re + im * im) / denom; // Re
            fourie_.in_[i][1] = (im * re - re * im) / denom; // Im
        }
        fourie_.transform(FFTW_FORWARD);

        double k = 2 * M_PI / Lambda_;
        for (int i = 0; i < height_; ++i) {
            double y = i - height_ / 2.0;
            for (int j = 0; j < width_; ++j) {
                double x = j - width_ / 2.0;
                // Фазовая модуляция
                double phase = copysign(pow(std::fabs(x), 2.5), x) + 
                    copysign(pow(std::fabs(y), 2.5), y);

                fourie_.in_[i * width_ + j][0] = fourie_.out_[i * width_ + j][0] / (cos(k * phase) + epsilon_);
                fourie_.in_[i * width_ + j][1] = fourie_.out_[i * width_ + j][1] / (sin(k * phase) + epsilon_) - fourie_.in_[i * width_ + j][0];
            }
        }
        fourie_.transform(FFTW_BACKWARD);

        return std::move(fourie_);
    }
};

// ==================================================
// Класс для отображения изображения через SFML
// ==================================================
class ImageRenderer {
private:
    sf::RenderWindow window_;   // Окно SFML
    sf::Texture texture_;       // Текстура для изображения
    sf::Sprite sprite_;         // Спрайт для отрисовки
    int width_;                 // Ширина изображения
    int height_;                // Высота изображения

public:
    ImageRenderer(int windowWidth, int windowHeight, int imgWidth, int imgHeight) 
        : window_(sf::VideoMode(windowWidth, windowHeight), "Imaging Simulation"),
          width_(imgWidth), height_(imgHeight) {
        texture_.create(imgWidth, imgHeight);
        sprite_.setTexture(texture_);
        // Масштабирование спрайта под размер окна
        sprite_.setScale(
            windowWidth / static_cast<float>(imgWidth), 
            windowHeight / static_cast<float>(imgHeight)
        );
    }

    // Обновление текстуры на основе данных интенсивности
    void updateImage(const fftw_complex *data) {

        sf::Uint8* pixels = new sf::Uint8[width_ * height_ * 4];
        // Find maximum value for normalization
        double max_val = 0.0;
        for (int i = 0; i < height_ * width_; ++i) {
            double val = sqrt(data[i][0] * data[i][0] + data[i][1] * data[i][1]);
            if (val > max_val) max_val = val;
        }

        // Convert to logarithmic scale for better visualization
        for (int i = 0; i < height_; ++i) {
            for (int j = 0; j < width_; ++j) {
                double val = sqrt(data[i * width_ + j][0] * data[i * width_ + j][0] + 
                            data[i * width_ + j][1] * data[i * width_ + j][1]);
                double log_val = log(1.0 + 100.0 * val / max_val);
                sf::Uint8 intensity = static_cast<sf::Uint8>(255 * log_val / log(101.0));

                pixels[(i * width_ + j) * 4] = intensity;
                pixels[(i * width_ + j) * 4 + 1] = intensity;
                pixels[(i * width_ + j) * 4 + 2] = intensity;
                pixels[(i * width_ + j) * 4 + 3] = 255;
            }
        }

        texture_.update(pixels);
        delete[] pixels;
    }

    // Отрисовка изображения
    void render() {
        window_.clear();
        window_.draw(sprite_);
        window_.display();
    }

    // Проверка открытости окна
    bool isOpen() const { return window_.isOpen(); }

    // Обработка событий
    bool pollEvent(sf::Event& event) { return window_.pollEvent(event); }

    // Закрытие окна
    void close() { window_.close(); }
};

// ==================================================
// Главная функция
// ==================================================
int main() {
    const int IMG_WIDTH = 1024;
    const int IMG_HEIGHT = 1024;
    const int LEN_X = 20;
    const int LEN_Y = 10;
    const double stepR = 0.5;
    const double stepLambda = 0.05e-6;
    double R = 15;
    double lambda = 0.5e-6;
    // CircleAperture system(IMG_WIDTH, IMG_HEIGHT, R, lambda);
    SerfaceEncoding system(IMG_WIDTH, IMG_HEIGHT, R, lambda);
    InverseFiltering system_restored(IMG_WIDTH, IMG_HEIGHT, lambda);
    ImageRenderer renderer(IMG_WIDTH, IMG_HEIGHT, IMG_WIDTH, IMG_HEIGHT);

    PointSource source(IMG_WIDTH, IMG_HEIGHT);
    // LineSource source(LEN_X);
    // CrossSource source(LEN_X, LEN_Y);

    const fftw_complex *intensity = system.computeIntencity(source);

    FourieTransform &&restored_fourie = system_restored.filter(intensity, system, source);

    CircleAperture new_system(restored_fourie, R, lambda);

    const fftw_complex *new_intensity = new_system.computeIntencity(source);
    renderer.updateImage(new_intensity);

    while (renderer.isOpen()) {
        sf::Event event;
        while (renderer.pollEvent(event)) {
            if (event.type == sf::Event::Closed) renderer.close();
            if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::J) new_system.changeR(stepR);
                else if (event.key.code == sf::Keyboard::K) new_system.changeR(-stepR);
                else if (event.key.code == sf::Keyboard::I) new_system.changeLambda(stepLambda);
                else if (event.key.code == sf::Keyboard::M) new_system.changeLambda(-stepLambda);
                new_intensity = new_system.computeIntencity(source);
                renderer.updateImage(new_intensity);
            }
        }

        renderer.render();
    }

    return 0;
}
