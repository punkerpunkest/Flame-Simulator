#pragma once
#include <SFML/Graphics.hpp>
#include <vector>

using GridData = std::vector<float>;

enum class BoundaryType{
    Scalar = 0,
    VelocityX = 1,
    VelocityY = 2,
};

constexpr int WINDOW_WIDTH = 800;
constexpr int WINDOW_HEIGHT = 600;

constexpr float DEFAULT_DIFFUSION = 0.001f;
constexpr float DEFAULT_VISCOSITY = 0.0001f;
constexpr int SOLVER_ITERATIONS = 20;

constexpr int computeGridSize(int width, int height, int resolution) noexcept {
    return ((width / resolution) + 2) * ((height / resolution) + 2);
}

class FluidSim {
public:
    explicit FluidSim(int res);
    ~FluidSim() = default;

    // fluid steps
    void velocityStep(float dt) noexcept;
    void densityStep(float dt) noexcept;

    // mouse input helpers
    void addDensityAt(float x, float y, float amount) noexcept;
    void addVelocityAt(float x, float y, float amountX, float amountY) noexcept;

    // clearing previous inputs
    void clearPrevDensity() noexcept;
    void clearPrevVelocity() noexcept;

    // drawing
    void displayDensity(sf::RenderWindow& window) const;
    void displayVelocity(sf::RenderWindow& window) const;

    // decay
    void decayDensity(float factor) noexcept;

private:
    int resolution_;
    int rows_;
    int N_;

    GridData dens_;
    GridData densPrev_;
    GridData u_;
    GridData v_;
    GridData uPrev_;
    GridData vPrev_;

    [[nodiscard]] int ix(int i, int j) const noexcept;

    void setBoundary(BoundaryType boundaryType, GridData& x) noexcept;
    void addSource(GridData& x, const GridData& s, float dt) noexcept;
    void diffuse(BoundaryType boundaryType, GridData& x, const GridData& x0,
                 float diffusionRate, float dt) noexcept;
    void advect(BoundaryType boundaryType, GridData& d, const GridData& d0,
                const GridData& velocityU, const GridData& velocityV, float dt) noexcept;
    void project(GridData& velocityU, GridData& velocityV,
                 GridData& pressure, GridData& divergence) noexcept;
};
