#include "FlameSim.hpp"
#include <algorithm>
#include <cmath>
#include <utility>
#include <SFML/Graphics.hpp>
FlameSim::FlameSim(int res)
    : resolution_{res}
    , rows_{WINDOW_HEIGHT / res}
    , N_{WINDOW_WIDTH / res}
    , dens_(computeGridSize(WINDOW_WIDTH, WINDOW_HEIGHT, res), 0.0f)
    , densPrev_(computeGridSize(WINDOW_WIDTH, WINDOW_HEIGHT, res), 0.0f)
    , u_(computeGridSize(WINDOW_WIDTH, WINDOW_HEIGHT, res), 0.0f)
    , v_(computeGridSize(WINDOW_WIDTH, WINDOW_HEIGHT, res), 0.0f)
    , uPrev_(computeGridSize(WINDOW_WIDTH, WINDOW_HEIGHT, res), 0.0f)
    , vPrev_(computeGridSize(WINDOW_WIDTH, WINDOW_HEIGHT, res), 0.0f)
    , temp_(computeGridSize(WINDOW_WIDTH, WINDOW_HEIGHT, res), AMBIENT_TEMPERATURE)
    , tempPrev_(computeGridSize(WINDOW_WIDTH, WINDOW_HEIGHT, res), 0.0f)
    , vorticity_(computeGridSize(WINDOW_WIDTH, WINDOW_HEIGHT, res), 0.0f)
    
{
}

int FlameSim::ix(int i, int j) const noexcept {
    return i + (N_ + 2) * j;
}

void FlameSim::setBoundary(BoundaryType boundaryType, GridData& x) noexcept {
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

void FlameSim::addSource(GridData& x, const GridData& s, float dt) noexcept {
    for (std::size_t i = 0; i < x.size(); ++i) {
        x[i] += dt * s[i];
    }
}

void FlameSim::diffuse(BoundaryType boundaryType, GridData& x, const GridData& x0, float diffusionRate, float dt) noexcept {
    const auto a = dt * diffusionRate * static_cast<float>(N_ * rows_);
    const auto divisor = 1.0f + 4.0f * a;

    for (int k = 0; k < SOLVER_ITERATIONS; ++k) {
        for (int i = 1; i <= N_; ++i) {
            for (int j = 1; j <= rows_; ++j) {
                x[ix(i,j)] = (x0[ix(i,j)] + a*(x[ix(i-1,j)] + x[ix(i+1,j)] + x[ix(i,j-1)] + x[ix(i,j+1)])) / divisor;
            }
        }
        setBoundary(boundaryType, x);
    }
}

void FlameSim::advect(BoundaryType boundaryType, GridData& d, const GridData& d0,
                      const GridData& velocityU, const GridData& velocityV, float dt) noexcept {
    const auto dt0 = dt * static_cast<float>(N_);
    const auto maxX = static_cast<float>(N_) + 0.5f;
    const auto maxY = static_cast<float>(rows_) + 0.5f;

    for (int i = 1; i <= N_; ++i) {
        for (int j = 1; j <= rows_; ++j) {
            auto x = static_cast<float>(i) - dt0 * velocityU[ix(i,j)];
            auto y = static_cast<float>(j) - dt0 * velocityV[ix(i,j)];

            x = std::max(0.5f, std::min(x, maxX));
            y = std::max(0.5f, std::min(y, maxY));

            const int i0 = static_cast<int>(x);
            const int j0 = static_cast<int>(y);
            const int i1 = i0 + 1;
            const int j1 = j0 + 1;

            const float s1 = x - static_cast<float>(i0);
            const float t1 = y - static_cast<float>(j0);
            const float s0 = 1.0f - s1;
            const float t0 = 1.0f - t1;

            d[ix(i,j)] = s0*(t0*d0[ix(i0,j0)] + t1*d0[ix(i0,j1)]) + s1*(t0*d0[ix(i1,j0)] + t1*d0[ix(i1,j1)]);
        }
    }
    setBoundary(boundaryType, d);
}

void FlameSim::project(GridData& velocityU, GridData& velocityV, GridData& pressure, GridData& divergence) noexcept {
    const auto h = 1.0f / static_cast<float>(N_);
    const auto halfH = 0.5f * h;

    for (int i=1; i<=N_; ++i){
        for (int j=1; j<=rows_; ++j){
            divergence[ix(i,j)] = -halfH*(velocityU[ix(i+1,j)] - velocityU[ix(i-1,j)] +
                                          velocityV[ix(i,j+1)] - velocityV[ix(i,j-1)]);
            pressure[ix(i,j)] = 0.0f;
        }
    }

    setBoundary(BoundaryType::Scalar, divergence);
    setBoundary(BoundaryType::Scalar, pressure);

    for(int k=0;k<SOLVER_ITERATIONS;++k){
        for(int i=1;i<=N_;++i){
            for(int j=1;j<=rows_;++j){
                pressure[ix(i,j)] = (divergence[ix(i,j)] + pressure[ix(i-1,j)] + pressure[ix(i+1,j)] + pressure[ix(i,j-1)] + pressure[ix(i,j+1)]) * 0.25f;
            }
        }
        setBoundary(BoundaryType::Scalar, pressure);
    }

    const auto halfOverH = 0.5f / h;
    for(int i=1;i<=N_;++i){
        for(int j=1;j<=rows_;++j){
            velocityU[ix(i,j)] -= halfOverH*(pressure[ix(i+1,j)] - pressure[ix(i-1,j)]);
            velocityV[ix(i,j)] -= halfOverH*(pressure[ix(i,j+1)] - pressure[ix(i,j-1)]);
        }
    }

    setBoundary(BoundaryType::VelocityX, velocityU);
    setBoundary(BoundaryType::VelocityY, velocityV);
}

void FlameSim::applyBuoyancy(float dt) noexcept {
    for (int i = 1; i <= N_; ++i) {
        for (int j = 1; j <= rows_; ++j) {
            const float density = dens_[ix(i,j)];
            const float temperature = temp_[ix(i,j)];

            const float buoyancyForce = BUOYANCY_SMOKE_WEIGHT * density 
                                      - BUOYANCY_HEAT_LIFT * (temperature - AMBIENT_TEMPERATURE);

            vPrev_[ix(i,j)] += dt * buoyancyForce;
        }
    }
}

void FlameSim::applyCooling(float dt) noexcept {
    for (int i = 1; i <= N_; ++i) {
        for (int j = 1; j <= rows_; ++j) {
            float& t = temp_[ix(i,j)];
            if (t > AMBIENT_TEMPERATURE) {
                const float deltaT = std::min(t - AMBIENT_TEMPERATURE, 5.0f);
                const float coolingRate = COOLING_COEFFICIENT * deltaT * deltaT * deltaT * deltaT;
                
                t -= coolingRate * dt;
                t = std::max(t, AMBIENT_TEMPERATURE);
            }
        }
    }
}

void FlameSim::applyVorticityConfinement(float dt) noexcept {
    for (int i = 1; i <= N_; ++i) {
        for (int j = 1; j <= rows_; ++j) {
            const float dvdx = (v_[ix(i+1,j)] - v_[ix(i-1,j)]) * 0.5f;
            const float dudy = (u_[ix(i,j+1)] - u_[ix(i,j-1)]) * 0.5f;
            vorticity_[ix(i,j)] = dvdx - dudy;
        }
    }
    
    for (int i = 2; i < N_; ++i) {
        for (int j = 2; j < rows_; ++j) {
            const float dωdx = (std::abs(vorticity_[ix(i+1,j)]) - std::abs(vorticity_[ix(i-1,j)])) * 0.5f;
            const float dωdy = (std::abs(vorticity_[ix(i,j+1)]) - std::abs(vorticity_[ix(i,j-1)])) * 0.5f;

            const float len = std::sqrt(dωdx * dωdx + dωdy * dωdy) + 1e-5f;
            const float nx = dωdx / len;
            const float ny = dωdy / len;

            const float ω = vorticity_[ix(i,j)];
            uPrev_[ix(i,j)] += dt * VORTICITY_EPSILON * ny * ω;
            vPrev_[ix(i,j)] -= dt * VORTICITY_EPSILON * nx * ω;
        }
    }
}


void FlameSim::velocityStep(float dt) noexcept {
    applyBuoyancy(dt);
    applyVorticityConfinement(dt);  

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

void FlameSim::densityStep(float dt) noexcept {
    addSource(dens_, densPrev_, dt);
    std::swap(dens_, densPrev_);
    diffuse(BoundaryType::Scalar, dens_, densPrev_, DEFAULT_DIFFUSION, dt);
    std::swap(dens_, densPrev_);
    advect(BoundaryType::Scalar, dens_, densPrev_, u_, v_, dt);
}

void FlameSim::temperatureStep(float dt) noexcept {
    addSource(temp_, tempPrev_, dt);
    
    std::swap(temp_, tempPrev_);
    diffuse(BoundaryType::Scalar, temp_, tempPrev_, DEFAULT_DIFFUSION, dt);
    
    std::swap(temp_, tempPrev_);
    advect(BoundaryType::Scalar, temp_, tempPrev_, u_, v_, dt);
    
    applyCooling(dt);
}

void FlameSim::clearPrevDensity() noexcept { std::fill(densPrev_.begin(), densPrev_.end(), 0.0f); }
void FlameSim::clearPrevVelocity() noexcept { std::fill(uPrev_.begin(), uPrev_.end(), 0.0f); std::fill(vPrev_.begin(), vPrev_.end(), 0.0f); }
void FlameSim::clearPrevTemperature() noexcept {std::fill(tempPrev_.begin(), tempPrev_.end(), 0.0f); }

void FlameSim::addDensityAt(float x, float y, float amount) noexcept {
    int i = static_cast<int>(x / resolution_) + 1;
    int j = static_cast<int>(y / resolution_) + 1;
    i = std::max(1, std::min(i, N_));
    j = std::max(1, std::min(j, rows_));
    densPrev_[ix(i,j)] = amount;
}

void FlameSim::addVelocityAt(float x, float y, float amountX, float amountY) noexcept {
    const int ci = static_cast<int>(x / resolution_) + 1;
    const int cj = static_cast<int>(y / resolution_) + 1;

    for (int di=-1; di<=1; ++di){
        for(int dj=-1;dj<=1;++dj){
            const int i = std::max(1,std::min(ci+di,N_));
            const int j = std::max(1,std::min(cj+dj,rows_));
            uPrev_[ix(i,j)] = amountX;
            vPrev_[ix(i,j)] = amountY;
        }
    }
}

void FlameSim::addTemperatureAt(float x, float y, float amount) noexcept {
    const int ci = static_cast<int>(x / resolution_) + 1;
    const int cj = static_cast<int>(y / resolution_) + 1;

    for (int di = -1; di <= 1; ++di) {
        for (int dj = -1; dj <= 1; ++dj) {
            const int i = std::max(1, std::min(ci + di, N_));
            const int j = std::max(1, std::min(cj + dj, rows_));
            tempPrev_[ix(i,j)] = amount;
        }
    }
}

void FlameSim::displayDensity(sf::RenderWindow& window) const {
    for(int i=1;i<=N_;++i){
        for(int j=1;j<=rows_;++j){
            float d = dens_[ix(i,j)];
            if(d > 0.001f){
                d = std::min(d*3.0f, 1.0f);
                sf::RectangleShape rect{{static_cast<float>(resolution_), static_cast<float>(resolution_)}};
                rect.setPosition(sf::Vector2f{
                    static_cast<float>((i-1)*resolution_),
                    static_cast<float>((j-1)*resolution_)
                });
                rect.setFillColor(sf::Color(255,255,255, static_cast<unsigned char>(d*255)));
                window.draw(rect);
            }
        }
    }
}

void FlameSim::displayVelocity(sf::RenderWindow& window) const {
    for(int i=1;i<=N_;++i){
        for(int j=1;j<=rows_;++j){
            const auto x = static_cast<float>((i-1)*resolution_);
            const auto y = static_cast<float>((j-1)*resolution_);
            const auto vx = u_[ix(i,j)];
            const auto vy = v_[ix(i,j)];
            const auto len = static_cast<float>(resolution_-2);

            sf::Vertex line[] = {
                {{x,y}, sf::Color(255,0,0,100)},
                {{x+vx*len, y+vy*len}, sf::Color(255,0,0,100)}
            };
            window.draw(line, 2, sf::PrimitiveType::Lines);
        }
    }
}

void FlameSim::displayFire(sf::RenderWindow& window) const {
     for (int i = 1; i <= N_; ++i) {
        for (int j = 1; j <= rows_; ++j) {
            const float d = dens_[ix(i,j)];
            const float t = temp_[ix(i,j)];
            
            if (d < 0.01f || t > 0.8f) continue;
            
            const float smokeAlpha = std::min(d * 0.5f, 0.4f) * (1.0f - t / 0.8f);
            const unsigned char gray = 40;

            sf::RectangleShape rect{{static_cast<float>(resolution_), static_cast<float>(resolution_)}};
            rect.setPosition(sf::Vector2f{
                static_cast<float>((i-1)*resolution_),
                static_cast<float>((j-1)*resolution_)
            });
            rect.setFillColor(sf::Color(gray, gray, gray, static_cast<unsigned char>(smokeAlpha * 255)));
            window.draw(rect);
        }
    }

    for (int i = 1; i <= N_; ++i) {
        for (int j = 1; j <= rows_; ++j) {
            const float t = temp_[ix(i,j)];
            if (t < 0.1f) continue;

            const float normalized = std::min(t / 2.0f, 1.0f);
            const float glowSize = resolution_ * (1.5f + normalized);
            const float r = std::min(1.0f, normalized * 1.5f);
            const float g = std::min(1.0f, normalized * 0.4f);
            const float b = 0.0f;
            const float a = normalized * 0.15f;

            sf::RectangleShape rect{{glowSize, glowSize}};
            rect.setOrigin(sf::Vector2f{glowSize * 0.25f, glowSize * 0.25f});
            rect.setPosition(sf::Vector2f{
                static_cast<float>((i-1)*resolution_),
                static_cast<float>((j-1)*resolution_)
            });
            rect.setFillColor(sf::Color(
                static_cast<unsigned char>(r * 255),
                static_cast<unsigned char>(g * 255),
                static_cast<unsigned char>(b * 255),
                static_cast<unsigned char>(a * 255)
            ));
            window.draw(rect, sf::BlendAdd);
        }
    }

    for (int i = 1; i <= N_; ++i) {
        for (int j = 1; j <= rows_; ++j) {
            const float t = temp_[ix(i,j)];
            if (t < 0.01f) continue;

            const float normalized = std::min(t / 2.0f, 1.0f);
            const float r = std::min(1.0f, normalized * 2.0f);
            const float g = std::min(1.0f, std::max(0.0f, normalized * normalized * 1.5f));
            const float b = std::min(1.0f, std::max(0.0f, (normalized - 0.5f) * normalized));
            const float a = std::min(1.0f, normalized * 0.9f);

            sf::RectangleShape rect{{static_cast<float>(resolution_), static_cast<float>(resolution_)}};
            rect.setPosition(sf::Vector2f{
                static_cast<float>((i-1)*resolution_),
                static_cast<float>((j-1)*resolution_)
            });
            rect.setFillColor(sf::Color(
                static_cast<unsigned char>(r * 255),
                static_cast<unsigned char>(g * 255),
                static_cast<unsigned char>(b * 255),
                static_cast<unsigned char>(a * 255)
            ));
            window.draw(rect, sf::BlendAdd);
        }
    }

    
    for (int i = 1; i <= N_; ++i) {
        for (int j = 1; j <= rows_; ++j) {
            const float t = temp_[ix(i,j)];
            if (t < 1.5f) continue;  

            const float intensity = std::min((t - 1.5f) / 2.0f, 1.0f);
            const float a = intensity * 0.6f;

            sf::RectangleShape rect{{static_cast<float>(resolution_), static_cast<float>(resolution_)}};
            rect.setPosition(sf::Vector2f{
                static_cast<float>((i-1)*resolution_),
                static_cast<float>((j-1)*resolution_)
            });
            rect.setFillColor(sf::Color(
                255,
                static_cast<unsigned char>(255 * intensity),
                static_cast<unsigned char>(200 * intensity),
                static_cast<unsigned char>(a * 255)
            ));
            window.draw(rect, sf::BlendAdd);
        }
    }
}

void FlameSim::decayDensity(float factor) noexcept {
    for(auto& d : dens_) d *= factor;
}
