#include <fftw3.h>
#include <SFML/Graphics.hpp>
#include <cmath>
#include <iostream>
#include <sstream>
#include "image.h"

int main() {
    const int N = 1024;
    double R = 10.0;
    const double R_step = 0.5;

    DotObject dot(N, N, R);


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

    dot.updateImage(R,texture);

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
                    dot.updateImage(R,texture);
                }
                else if (event.key.code == sf::Keyboard::K) {
                    R = std::max(R - R_step, R_step);
                    dot.updateImage(R,texture);
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

    return 0;
}
