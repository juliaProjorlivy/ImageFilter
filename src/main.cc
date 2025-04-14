#include <fftw3.h>
#include <SFML/Graphics.hpp>
#include <cmath>
#include <iostream>
#include <sstream>

void updateImage(int N, double R, sf::Uint8* pixels, sf::Texture& texture) {
    // Массив для функции пропускания
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);

    // Заполнение функции пропускания (круглая отверстие -> дельта-функция Дирака)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double x = i - N / 2.0;
            double y = j - N / 2.0;
            double r = sqrt(x * x + y * y);
            in[i * N + j][0] = (r <= R) ? 1.0 : 0.0;
            in[i * N + j][1] = 0.0;
        }
    }

    fftw_plan plan = fftw_plan_dft_2d(N, N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Calculate I
    double max_intensity = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double real = out[i * N + j][0];
            double imag = out[i * N + j][1];
            double intensity = real * real + imag * imag;
            if (intensity > max_intensity) max_intensity = intensity;
        }
    }

    // Normilization and image filling
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double real = out[i * N + j][0];
            double imag = out[i * N + j][1];
            double intensity = (real * real + imag * imag) / max_intensity * 255;
            sf::Uint8 value = static_cast<sf::Uint8>(intensity);
            pixels[(i * N + j) * 4] = value;     // R
            pixels[(i * N + j) * 4 + 1] = value; // G
            pixels[(i * N + j) * 4 + 2] = value; // B
            pixels[(i * N + j) * 4 + 3] = 255;   // A
        }
    }

    texture.update(pixels);

    // Free everything
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
}

int main() {
    const int N = 512;
    double R = 10.0;
    const double R_step = 0.5;

    sf::Uint8 *pixels = new sf::Uint8[N * N * 4]; // RGBA array

    sf::Texture texture;
    texture.create(N, N);
    sf::Sprite sprite(texture);

    sf::Font font;
    if (!font.loadFromFile("../arial.ttf")) {
        std::cerr << "Ошибка загрузки шрифта!" << std::endl;
        return 1;
    }

    sf::Text radiusText;
    radiusText.setFont(font);
    radiusText.setCharacterSize(20);
    radiusText.setFillColor(sf::Color::Red);
    radiusText.setPosition(10, 10);

    updateImage(N, R, pixels, texture);

    sf::RenderWindow window(sf::VideoMode(N, N), "Airy diffraction pattern");
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();

            // j - to increase R, k - to decrease R
            if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::J) {
                    R += R_step;
                    updateImage(N, R, pixels, texture);
                }
                else if (event.key.code == sf::Keyboard::K) {
                    R = std::max(R - R_step, R_step);
                    updateImage(N, R, pixels, texture);
                }
            }
        }

        std::ostringstream oss;
        oss << R;
        radiusText.setString(oss.str());

        window.clear();
        window.draw(sprite);
        window.draw(radiusText);
        window.display();
    }

    delete[] pixels;
    return 0;
}
