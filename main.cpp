#include <SFML/Graphics.hpp>
#include "FlameSim.hpp"
#include <ctime>
#include <cstdlib>
#include <iostream>

int main() {

    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    sf::RenderWindow window{sf::VideoMode{{WINDOW_WIDTH, WINDOW_HEIGHT}}, "Flame Simulator"};
    window.setFramerateLimit(60);

    FlameSim flame{10};

    constexpr float dt = 0.05f;
    constexpr float velocityScale = 2.0f;
    constexpr float densityAmount = 5.0f;
    constexpr float temperatureAmount = 5.0f;
    constexpr float decayFactor = 0.995f;
    

    sf::Vector2i lastMouse{0, 0};
    bool firstMouse = true;

    sf::Clock printClock;

    while(window.isOpen()){
        while(auto event = window.pollEvent()){
            if(event->is<sf::Event::Closed>()){
                window.close();
            }
        }

        flame.clearPrevDensity();
        flame.clearPrevVelocity();
        flame.clearPrevTemperature();

        if(sf::Mouse::isButtonPressed(sf::Mouse::Button::Left)){
            const auto mousePos = sf::Mouse::getPosition(window);
            if(!firstMouse){
                const auto dx = static_cast<float>(mousePos.x - lastMouse.x);
                const auto dy = static_cast<float>(mousePos.y - lastMouse.y);
                flame.addVelocityAt(
                    static_cast<float>(mousePos.x),
                    static_cast<float>(mousePos.y),
                    dx*velocityScale,
                    dy*velocityScale
                );
            }
            flame.addDensityAt(
                static_cast<float>(mousePos.x),
                static_cast<float>(mousePos.y),
                densityAmount
            );
            flame.addTemperatureAt(
                static_cast<float>(mousePos.x),
                static_cast<float>(mousePos.y),
                temperatureAmount  
            );

            lastMouse = mousePos;
            firstMouse = false;
        } else {
            firstMouse = true;
        }

        flame.velocityStep(dt);
        flame.densityStep(dt);
        flame.temperatureStep(dt);
        flame.decayDensity(decayFactor);

        window.clear(sf::Color::Black);
        flame.displayFire(window);
        window.display();

        if (printClock.getElapsedTime().asSeconds() >= 5.0f) {
            Profiler::get().print();
            printClock.restart();
        }
    }

    return 0;
}
