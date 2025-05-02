#include "../include/ImageRenderer.h"
#include "../include/OpticalSystems.h"
#include "../include/LightSources.h"
#include "../include/InverseFiltering.h"
#include <SFML/Window/Keyboard.hpp>
#include <stdio.h>

int main()
{
    const int IMG_SIZE = 512;
    const int WINDOW_SIZE = 1024;
    const int LEN_X = 30;
    const int LEN_Y = 50;
    double R = 10;
    double lambda = 0.5e-6;
    const double stepR = 0.5;
    const double stepLambda = 0.05e-6;
    const int stepX = 10;
    const int stepY = 10;
    
    // Инициализация компонентов
    ImageRenderer renderer(WINDOW_SIZE, IMG_SIZE);
    
    // Создание оптической системы
    CircleAperture circleAperture(IMG_SIZE, IMG_SIZE, R, lambda);
    SerfaceEncoding destorer(IMG_SIZE, IMG_SIZE, R, lambda);
    InverseFiltering restorer(IMG_SIZE, IMG_SIZE, lambda);
    
    // Источники света

    int counter = 0;
    PointSource pointSource(IMG_SIZE, IMG_SIZE);
    LineSource lineSource(IMG_SIZE, IMG_SIZE, LEN_X);
    CrossSource crossSource(IMG_SIZE, IMG_SIZE, LEN_X, LEN_Y);

    std::vector<LightSource *> sources = {&pointSource, &lineSource, &crossSource};
    LightSource *source = sources[counter];

    // Основной цикл
    while (renderer.isOpen())
    {
        sf::Event event;
        while (renderer.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                renderer.close();
            
            if (event.type == sf::Event::KeyPressed)
            {
                switch(event.key.code) {
                    case sf::Keyboard::J: 
                    {
                        circleAperture.changeR(stepR);
                        destorer.changeR(stepR);
                        break;
                    }
                    case sf::Keyboard::K:
                    {
                        circleAperture.changeR(-stepR);
                        destorer.changeR(-stepR);
                        break;
                    }
                    case sf::Keyboard::I:
                    {
                        circleAperture.changeLambda(stepLambda);
                        destorer.changeLambda(stepLambda);
                        break;
                    }
                    case sf::Keyboard::M:
                    {
                        circleAperture.changeLambda(-stepLambda);
                        destorer.changeLambda(-stepLambda);
                        break;
                    }
                    case sf::Keyboard::Equal:
                    {
                        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Y))
                        {
                            if (counter % sources.size() == 2)
                            {
                                source->changeLenY(stepY);
                            }
                        }
                        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::X))
                        {
                            if (counter % sources.size() != 0)
                            {
                                source->changeLenX(stepX);
                            }
                        }

                        break;
                    }
                    case sf::Keyboard::Hyphen:
                    {
                        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Y))
                        {
                            if (counter % sources.size() == 2)
                            {
                                source->changeLenY(-stepY);
                            }
                        }
                        else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Hyphen))
                        {
                            if (counter % sources.size() != 0)
                            {
                                source->changeLenX(-stepX);
                            }
                        }

                        break;
                    }
                    case sf::Keyboard::Escape:
                    {
                        renderer.close();
                        break;
                    }
                    case sf::Keyboard::Enter:
                    {
                        source = sources[(++counter) % sources.size()];
                        break;
                    }
                    case sf::Keyboard::Q:
                    {
                        source = &pointSource;
                        counter = 0;
                        break;
                    }
                    case sf::Keyboard::W:
                    {
                        source = &lineSource;
                        counter = 1;
                        break;
                    }
                    case sf::Keyboard::E:
                    {
                        source = &crossSource;
                        counter = 2;
                        break;
                    }
                    default: 
                        break;
                }
            }
        }

        // Получение данных для отображения
        const fftw_complex *intensity = circleAperture.computeIntencity(*source);
        const fftw_complex *dest_intensity = destorer.computeIntencity(*source);
        
        FourieTransform &&restored_fourie = restorer.filter(intensity, circleAperture, *source);
        CircleAperture new_system(restored_fourie, circleAperture.getR(), circleAperture.getLambda());
        const fftw_complex *new_intensity = new_system.computeIntencity(*source);
        
        // Обновление отображения
        renderer.updateImage(0, intensity);     // Левый верхний - исходное
        renderer.updateImage(1, dest_intensity);    // Правый верхний - искаженное
        renderer.updateImage(3, new_intensity);// Правый нижний - восстановленное
        renderer.updateInfo(circleAperture.getR(), circleAperture.getSize(), circleAperture.getLambda(), counter, source->getLenX(), source->getLenY());
        
        renderer.render();
    }

    return 0;
}