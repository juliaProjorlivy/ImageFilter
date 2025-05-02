#pragma once
#include <SFML/Graphics.hpp>
#include <string>
#include <fftw3.h>
#include <cmath>

class ImageRenderer
{
private:
    sf::RenderWindow window_;
    sf::Texture textures_[4];
    sf::Sprite sprites_[4];
    sf::Font font_;
    sf::Text infoText_;
    int renderSize_;
    int fullWidth_;
    int fullHeight_;

    std::vector<std::string> LightSourceName = {"Point source",
                                                "Line source",
                                                "Cross source"};

    void setupSprite(sf::Sprite& sprite, int index)
    {
        float scaleX = (fullWidth_ / 2.0f) / textures_[index].getSize().x;
        float scaleY = (fullHeight_ / 2.0f) / textures_[index].getSize().y;
        sprite.setScale(scaleX, scaleY);
        sprite.setPosition((index % 2) * (fullWidth_ / 2), 
                          (index / 2) * (fullHeight_ / 2));
    }

public:
    ImageRenderer(int windowSize, int renderSize) 
        : window_(sf::VideoMode(windowSize, windowSize), "Optics Simulation"),
          renderSize_(renderSize),
          fullWidth_(windowSize),
          fullHeight_(windowSize)
    {
        
        // Инициализация текстур
        for(auto& tex : textures_)
        {
            tex.create(renderSize_, renderSize_);
        }

        // Настройка спрайтов
        for(int i = 0; i < 4; ++i)
        {
            sprites_[i].setTexture(textures_[i]);
            setupSprite(sprites_[i], i);
        }

        // Настройка текста
        if(!font_.loadFromFile("arial.ttf"))
        {
            throw std::runtime_error("Failed to load font!");
        }

        infoText_.setFont(font_);
        infoText_.setCharacterSize(20);
        infoText_.setFillColor(sf::Color::White);
        infoText_.setPosition(20, fullHeight_ - 200);
    }

    void updateImage(int index, const fftw_complex* data)
    {
        sf::Uint8* pixels = new sf::Uint8[renderSize_ * renderSize_ * 4];
        
        // Поиск максимального значения для нормализации
        double maxVal = 0.0;
        
        for(int i = 0; i < renderSize_ * renderSize_; ++i)
        {
            double val = sqrt(data[i][0]*data[i][0] + data[i][1]*data[i][1]);
            if(val > maxVal) maxVal = val;
        }

        // Заполнение пикселей с логарифмической шкалой
        for(int i = 0; i < renderSize_ * renderSize_; ++i)
        { 
            double val = sqrt(data[i][0]*data[i][0] + data[i][1]*data[i][1]);
            double logVal = log(1.0 + 100.0 * val / maxVal);
            sf::Uint8 intensity = static_cast<sf::Uint8>(255 * logVal / log(101.0));
            
            pixels[i*4]   = intensity;
            pixels[i*4+1] = intensity;
            pixels[i*4+2] = intensity;
            pixels[i*4+3] = 255;
        }

        textures_[index].update(pixels);
        delete[] pixels;
    }

    void updateInfo(float R, int size, float lambda, int counter, int lenX, int lenY)
    {
        int index = counter % LightSourceName.size();

        std::string info = "Parameters:\n";
        info += "Radius: " + std::to_string(R) + " px\n";
        info += "Image Size: " + std::to_string(size) + "x" + std::to_string(size) + "\n";
        info += "Wavelength: " + std::to_string(lambda*1e9) + " nm\n";
        info += "LightSource: " + LightSourceName[index] + "\n";
        
        if (index != 0)
        {
            info += "Length X: " + std::to_string(lenX) + "\n";
        }

        if (index == 2)
        {
            info += "Length Y: " + std::to_string(lenY);
        }

        infoText_.setString(info);
    }

    void render()
    {
        window_.clear();

        for(auto& sprite : sprites_)
        {
            window_.draw(sprite);
        }

        window_.draw(infoText_);
        window_.display();
    }

    bool isOpen() const 
    { 
        return window_.isOpen();
    
    }
    bool pollEvent(sf::Event& event) 
    { 
        return window_.pollEvent(event); 
    }

    void close() 
    { 
        window_.close();
    }
};