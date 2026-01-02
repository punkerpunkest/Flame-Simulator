#include <SFML/Graphics.hpp>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>

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

class FluidSim{
public:
  explicit FluidSim(int res)
    : resolution_{res}
    , rows_{WINDOW_HEIGHT / res}
    , N_{WINDOW_WIDTH / res}
    , dens_(computeGridSize(WINDOW_WIDTH, WINDOW_HEIGHT, res), 0.0f)
    , densPrev_(computeGridSize(WINDOW_WIDTH, WINDOW_HEIGHT, res), 0.0f)
    , u_(computeGridSize(WINDOW_WIDTH, WINDOW_HEIGHT, res), 0.0f)
    , v_(computeGridSize(WINDOW_WIDTH, WINDOW_HEIGHT, res), 0.0f)
    , uPrev_(computeGridSize(WINDOW_WIDTH, WINDOW_HEIGHT, res), 0.0f)
    , vPrev_(computeGridSize(WINDOW_WIDTH, WINDOW_HEIGHT, res), 0.0f)
  {

  }

  ~FluidSim() = default;
  FluidSim(const FluidSim&) = default;
  FluidSim& operator=(const FluidSim&) = default;
  FluidSim(FluidSim&&) noexcept = default;
  FluidSim& operator=(FluidSim&&) noexcept = default;


  [[nodiscard]] int ix(int i, int j) const noexcept {
    return i + (N_ + 2) * j;
  }
  
  void setBoundary(BoundaryType boundaryType, GridData& x) noexcept{
    const bool reflectHorizontal = (boundaryType == BoundaryType::VelocityX);
    const bool reflectVertical = (boundaryType == BoundaryType::VelocityY);

    const float hFactor = reflectHorizontal ? -1.0f : 1.0f;
    const float vFactor = reflectVertical ? -1.0f : 1.0f;

    for (int i = 1; i <= N_; ++i) {
      
      x[ix(0, i)]      = hFactor * x[ix(1, i)];
      x[ix(N_ + 1, i)] = hFactor * x[ix(N_, i)];

      x[ix(i, 0)]         = vFactor * x[ix(i, 1)];
      x[ix(i, rows_ + 1)] = vFactor * x[ix(i, rows_)];
    }

    x[ix(0, 0)]             = 0.5f * (x[ix(1, 0)] + x[ix(0, 1)]);
    x[ix(0, rows_ + 1)]     = 0.5f * (x[ix(1, rows_ + 1)] + x[ix(0, rows_)]);
    x[ix(N_ + 1, 0)]        = 0.5f * (x[ix(N_, 0)] + x[ix(N_ + 1, 1)]);
    x[ix(N_ + 1, rows_ + 1)] = 0.5f * (x[ix(N_, rows_ + 1)] + x[ix(N_ + 1, rows_)]);
  }

  void addSource(GridData& x, const GridData& s, float dt) noexcept {
    const auto size = x.size();
    for (std::size_t i = 0; i < size; ++i) {
      x[i] += dt * s[i];
    }
  }

  void diffuse(BoundaryType boundaryType, GridData& x, const GridData& x0, 
               float diffusionRate, float dt) noexcept {

    const auto a = dt * diffusionRate * static_cast<float>(N_ * rows_);
    const auto divisor = 1.0f + 4.0f * a;
    
    for (int k = 0; k < SOLVER_ITERATIONS; ++k) {
      for (int i = 1; i <= N_; ++i) {
          for (int j = 1; j <= rows_; ++j) {
              x[ix(i, j)] = (x0[ix(i, j)] + a * (
                  x[ix(i - 1, j)] +
                  x[ix(i + 1, j)] +
                  x[ix(i, j - 1)] +
                  x[ix(i, j + 1)]
              )) / divisor;
          }
      }
      setBoundary(boundaryType, x);
    }
  }
    
  void advect(BoundaryType boundaryType, GridData& d, const GridData& d0,
              const GridData& velocityU, const GridData& velocityV, float dt) noexcept {
  
    const auto dt0 = dt * static_cast<float>(N_);
    const auto maxX = static_cast<float>(N_) + 0.5f;
    const auto maxY = static_cast<float>(rows_) + 0.5f;
  
    for (int i = 1; i <= N_; ++i) {
      for (int j = 1; j <= rows_; ++j) {
       
        auto x = static_cast<float>(i) - dt0 * velocityU[ix(i, j)];
        auto y = static_cast<float>(j) - dt0 * velocityV[ix(i, j)];

        x = std::max(0.5f, std::min(x, maxX));
        y = std::max(0.5f, std::min(y, maxY));

        const auto i0 = static_cast<int>(x);
        const auto j0 = static_cast<int>(y);
        const auto i1 = i0 + 1;
        const auto j1 = j0 + 1;

        const auto s1 = x - static_cast<float>(i0);
        const auto t1 = y - static_cast<float>(j0);
        const auto s0 = 1.0f - s1;
        const auto t0 = 1.0f - t1;

        d[ix(i, j)] = s0 * (t0 * d0[ix(i0, j0)] + t1 * d0[ix(i0, j1)]) +
                      s1 * (t0 * d0[ix(i1, j0)] + t1 * d0[ix(i1, j1)]);
      }
    }
    setBoundary(boundaryType, d);
  }

  void project(GridData& velocityU, GridData& velocityV,
              GridData& pressure, GridData& divergence) noexcept {
    
    const auto h = 1.0f / static_cast<float>(N_);
    const auto halfH = 0.5f * h;
      
    for (int i = 1; i <= N_; ++i) {
        for (int j = 1; j <= rows_; ++j) {
            divergence[ix(i, j)] = -halfH * (
                velocityU[ix(i + 1, j)] - velocityU[ix(i - 1, j)] +
                velocityV[ix(i, j + 1)] - velocityV[ix(i, j - 1)]
            );
            pressure[ix(i, j)] = 0.0f;
        }
    }

    setBoundary(BoundaryType::Scalar, divergence);
    setBoundary(BoundaryType::Scalar, pressure);

    for (int k = 0; k < SOLVER_ITERATIONS; ++k) {
      for (int i = 1; i <= N_; ++i) {
        for (int j = 1; j <= rows_; ++j) {
          pressure[ix(i, j)] = (
          divergence[ix(i, j)] +
          pressure[ix(i - 1, j)] +
          pressure[ix(i + 1, j)] +
          pressure[ix(i, j - 1)] +
          pressure[ix(i, j + 1)]
          ) * 0.25f;  
        }
      }
      setBoundary(BoundaryType::Scalar, pressure);
    }


    const auto halfOverH = 0.5f / h;
    for (int i = 1; i <= N_; ++i) {
      for (int j = 1; j <= rows_; ++j) {
        velocityU[ix(i, j)] -= halfOverH * (pressure[ix(i + 1, j)] - pressure[ix(i - 1, j)]);
        velocityV[ix(i, j)] -= halfOverH * (pressure[ix(i, j + 1)] - pressure[ix(i, j - 1)]);
      }
    }
    setBoundary(BoundaryType::VelocityX, velocityU);
    setBoundary(BoundaryType::VelocityY, velocityV);
  }

  void velocityStep(float dt) noexcept {
    addSource(u_, uPrev_, dt);
    addSource(v_, vPrev_, dt);

    std::swap(u_, uPrev_);
    diffuse(BoundaryType::VelocityX, u_, uPrev_, DEFAULT_VISCOSITY, dt);
    
    std::swap(v_, vPrev_);
    diffuse(BoundaryType::VelocityY, v_, vPrev_, DEFAULT_VISCOSITY, dt);
    
    project(u_, v_, uPrev_, vPrev_);
    
    std::swap(u_, uPrev_);
    std::swap(v_, vPrev_);
    
    advect(BoundaryType::VelocityX, u_, uPrev_, uPrev_, vPrev_, dt);
    advect(BoundaryType::VelocityY, v_, vPrev_, uPrev_, vPrev_, dt);
    
    project(u_, v_, uPrev_, vPrev_);
  }

  void densityStep(float dt) noexcept {
    addSource(dens_, densPrev_, dt);
    
    std::swap(dens_, densPrev_);
    diffuse(BoundaryType::Scalar, dens_, densPrev_, DEFAULT_DIFFUSION, dt);
    
    std::swap(dens_, densPrev_);
    advect(BoundaryType::Scalar, dens_, densPrev_, u_, v_, dt);
  }

  void clearPrevDensity() noexcept {
    std::fill(densPrev_.begin(), densPrev_.end(), 0.0f);
  }

  void clearPrevVelocity() noexcept {
    std::fill(uPrev_.begin(), uPrev_.end(), 0.0f);
    std::fill(vPrev_.begin(), vPrev_.end(), 0.0f);
  }

  void addDensityAt(float x, float y, float amount) noexcept {
    auto i = static_cast<int>(x / resolution_) + 1;
    auto j = static_cast<int>(y / resolution_) + 1;
    i = std::max(1, std::min(i, N_));
    j = std::max(1, std::min(j, rows_));
    densPrev_[ix(i, j)] = amount;
  }

  void addVelocityAt(float x, float y, float amountX, float amountY) noexcept {
    const auto ci = static_cast<int>(x / resolution_) + 1;
    const auto cj = static_cast<int>(y / resolution_) + 1;

    for (int di = -1; di <= 1; ++di) {
      for (int dj = -1; dj <= 1; ++dj) {
        const auto i = std::max(1, std::min(ci + di, N_));
        const auto j = std::max(1, std::min(cj + dj, rows_));
        uPrev_[ix(i, j)] = amountX;
        vPrev_[ix(i, j)] = amountY;
      }
    }
  }

  void displayVelocity(sf::RenderWindow& window) const {
    for (int i = 1; i <= N_; ++i) {
      for (int j = 1; j <= rows_; ++j) {
        const auto x = static_cast<float>((i - 1) * resolution_);
        const auto y = static_cast<float>((j - 1) * resolution_);
        const auto vx = u_[ix(i, j)];
        const auto vy = v_[ix(i, j)];
        const auto len = static_cast<float>(resolution_ - 2);
        
        sf::Vertex line[] = {
            {{x, y}, sf::Color(255, 0, 0, 100)},
            {{x + vx * len, y + vy * len}, sf::Color(255, 0, 0, 100)}
        };
        window.draw(line, 2, sf::PrimitiveType::Lines);
      }
    }
  }

  void displayDensity(sf::RenderWindow& window) const {
    for (int i = 1; i <= N_; ++i) {
      for (int j = 1; j <= rows_; ++j) {
        auto d = dens_[ix(i, j)];
        if (d > 0.001f) {
            d = std::min(d * 3.0f, 1.0f);
            sf::RectangleShape rect{{static_cast<float>(resolution_), 
                                     static_cast<float>(resolution_)}};
            rect.setPosition       ({static_cast<float>((i - 1) * resolution_), 
                                     static_cast<float>((j - 1) * resolution_)});
            rect.setFillColor       (sf::Color(255, 255, 255, 
                                     static_cast<unsigned char>(d * 255)));
            window.draw(rect);
        }
      }
    }
  }

  void decayDensity(float factor) noexcept {
    for (auto& d : dens_) {
        d *= factor;
    }
  }

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
};

int main() {
  std::srand(static_cast<unsigned int>(std::time(nullptr)));

  sf::RenderWindow window{
      sf::VideoMode{{WINDOW_WIDTH, WINDOW_HEIGHT}}, 
      "Fluid Simulator"
  };
  window.setFramerateLimit(60);

  FluidSim fluid{10};
  
  constexpr float dt = 0.05f;
  constexpr float velocityScale = 2.0f;
  constexpr float densityAmount = 5.0f;
  constexpr float decayFactor = 0.995f;


  sf::Vector2i lastMouse{0, 0};
  auto firstMouse = true;

  while (window.isOpen()){
    while(auto event = window.pollEvent()){
      if (event->is<sf::Event::Closed>()) {
        window.close();
      }
    }

    fluid.clearPrevDensity();
    fluid.clearPrevVelocity();

    if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Left)) {  
      const auto mousePos = sf::Mouse::getPosition(window);
      if (!firstMouse) {
        const auto dx = static_cast<float>(mousePos.x - lastMouse.x);
        const auto dy = static_cast<float>(mousePos.y - lastMouse.y);
        fluid.addVelocityAt(
          static_cast<float>(mousePos.x), 
          static_cast<float>(mousePos.y),
          dx * velocityScale, 
          dy * velocityScale
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