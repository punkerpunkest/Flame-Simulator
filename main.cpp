#include <SFML/Graphics.hpp>
#include "FluidSim.hpp"
#include <ctime>
#include <cstdlib>

int main() {
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    sf::RenderWindow window{sf::VideoMode{{WINDOW_WIDTH, WINDOW_HEIGHT}}, "Fluid Simulator"};
    window.setFramerateLimit(60);

    FluidSim fluid{10};

    constexpr float dt = 0.05f;
    constexpr float velocityScale = 2.0f;
    constexpr float densityAmount = 5.0f;
    constexpr float decayFactor = 0.995f;

    sf::Vector2i lastMouse{0, 0};
    bool firstMouse = true;

    while(window.isOpen()){
        while(auto event = window.pollEvent()){
            if(event->is<sf::Event::Closed>()){
                window.close();
            }
        }

        fluid.clearPrevDensity();
        fluid.clearPrevVelocity();

        if(sf::Mouse::isButtonPressed(sf::Mouse::Button::Left)){
            const auto mousePos = sf::Mouse::getPosition(window);
            if(!firstMouse){
                const auto dx = static_cast<float>(mousePos.x - lastMouse.x);
                const auto dy = static_cast<float>(mousePos.y - lastMouse.y);
                fluid.addVelocityAt(
                    static_cast<float>(mousePos.x),
                    static_cast<float>(mousePos.y),
                    dx*velocityScale,
                    dy*velocityScale
                );
            }
            fluid.addDensityAt(
                static_cast<float>(mousePos.x),
                static_cast<float>(mousePos.y),
                densityAmount
            );

            lastMouse = mousePos;
            firstMouse = false;
        } else {
            firstMouse = true;
        }

        fluid.velocityStep(dt);
        fluid.densityStep(dt);
        fluid.decayDensity(decayFactor);

        window.clear(sf::Color::Black);
        fluid.displayDensity(window);
        window.display();
    }

    return 0;
}
