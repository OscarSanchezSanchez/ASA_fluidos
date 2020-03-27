#include "Scene.h"

#include "Numeric/PCGSolver.h"

namespace asa
{
namespace
{
////////////////////////////////////////////////
// Add any reusable classes or functions HERE //
////////////////////////////////////////////////

    float checkInside(Index2 &a, Array2<float> &b)
    {
        Index2 sz = b.getSize();
        if (a.x >= 0 && a.x < sz.x && a.y >= 0 && a.y < sz.y)
            return b[a];
        else
            return 0.0f;
    }
}  // namespace

// advection
void Fluid2::fluidAdvection(const float dt)
{
    AABox2 box = grid.getDomain();
    {
        // Ink advecion HERE
        Array2<Vector3> ink(inkRGB);
        for (size_t i = 0; i < grid.getSize().x; i++)
        {
            for (size_t j = 0; j < grid.getSize().y; j++)
            {
                Index2 index_ij(i, j);
                Index2 pointX(i + 1, j);
                Index2 pointY(i, j + 1);

                float velX = 0.5f * (velocityX[index_ij] + velocityX[pointX]);
                float velY = 0.5f * (velocityY[index_ij] + velocityY[pointY]);
                Vector2 v(velX, velY);

                Vector2 prevPoint = grid.getCellPos(index_ij) - dt * v;

                if (prevPoint.x < box.minPosition.x)
                    prevPoint.x = box.minPosition.x;
                else if (prevPoint.x > box.maxPosition.x)
                    prevPoint.x = box.maxPosition.x;

                if (prevPoint.y < box.minPosition.y)
                    prevPoint.y = box.minPosition.y;
                else if (prevPoint.y > box.maxPosition.y)
                    prevPoint.y = box.maxPosition.y;

                Vector2 indexPrev = grid.getCellIndex(prevPoint);
                Vector2 min(floor(indexPrev.x), floor(indexPrev.y));
                Vector2 max(ceil(indexPrev.x), ceil(indexPrev.y));

                int minX = clamp(min.x, 0, inkRGB.getSize().x - 1);
                int maxX = clamp(max.x, 0, inkRGB.getSize().x - 1);
                int minY = clamp(min.y, 0, inkRGB.getSize().y - 1);
                int maxY = clamp(max.y, 0, inkRGB.getSize().y - 1);

                Index2 p1(minX, minY);
                Index2 p2(maxX, minY);
                Index2 p3(maxX, maxY);
                Index2 p4(minX, maxY);

                Vector2 diff = indexPrev - min;
                float alpha = diff.x;
                float betha = diff.y;
                Vector3 inkInterpolate = bilerp(ink[p1], ink[p2], ink[p3], ink[p4], alpha, betha);

                //aqui se deberia comprobar que no se salga
                inkRGB[index_ij] = inkInterpolate; 

            }
        }
    }

    {
        // Velocity acvection term HERE
        //if (Scene::testcase >= Scene::TEST_ADVECTION) 
        {
            Index2 faceSizeX = grid.getSizeFacesX();
            Index2 faceSizeY = grid.getSizeFacesY();
            Array2<float> VX_result(velocityX.getSize());
            Array2<float> VY_result(velocityY.getSize());


            for (int i = 0; i < faceSizeX.x; i++) 
            {
                for (int j = 0; j < faceSizeX.y; j++) 
                {
                    Index2 index_ij(i, j);
                    Vector2 point = grid.getFaceXPos(index_ij);

                    float vx = checkInside(index_ij, velocityX);

                    float vy0 = checkInside(index_ij, velocityY);
                    float vy1 = checkInside(Index2(index_ij.x, index_ij.y + 1), velocityY);
                    float vy2 = checkInside(Index2(index_ij.x - 1, index_ij.y + 1), velocityY);
                    float vy3 = checkInside(Index2(index_ij.x - 1, index_ij.y), velocityY);
                    float vy = bilerp(vy0, vy1, vy2, vy3, 0.5f, 0.5f);

                    Vector2 v(vx, vy);

                    Vector2 prevPoint = point - dt * v;

                    if (prevPoint.x < box.minPosition.x)
                        prevPoint.x = box.minPosition.x;
                    else if (prevPoint.x > box.maxPosition.x)
                        prevPoint.x = box.maxPosition.x;

                    if (prevPoint.y < box.minPosition.y)
                        prevPoint.y = box.minPosition.y;
                    else if (prevPoint.y > box.maxPosition.y)
                        prevPoint.y = box.maxPosition.y;

                    Vector2 indexF = grid.getFaceIndex(prevPoint, 0);
                    Vector2 min(floor(indexF.x), floor(indexF.y));
                    Vector2 max(ceil(indexF.x), ceil(indexF.y));

                    int minX = clamp(min.x, 0, inkRGB.getSize().x - 1);
                    int maxX = clamp(max.x, 0, inkRGB.getSize().x - 1);
                    int minY = clamp(min.y, 0, inkRGB.getSize().y - 1);
                    int maxY = clamp(max.y, 0, inkRGB.getSize().y - 1);

                    Index2 p1(minX, minY);
                    Index2 p2(maxX, minY);
                    Index2 p3(maxX, maxY);
                    Index2 p4(minX, maxY);

                    Vector2 diff = indexF - min;
                    float alpha = diff.x;
                    float betha = diff.y;
                    float vInterpolate = bilerp(velocityX[p1], velocityX[p2], velocityX[p4], velocityX[p3], alpha, betha);

                    VX_result[index_ij] = vInterpolate; 
                }
            }

            for (int i = 0; i < faceSizeY.x; i++) 
            {
                for (int j = 0; j < faceSizeY.y; j++) 
                {
                    Index2 index_ij(i, j);
                    Vector2 point = grid.getFaceYPos(index_ij);

                    float vy = checkInside(index_ij, velocityY);

                    float vx0 = checkInside(index_ij, velocityX);
                    float vx1 = checkInside(Index2(index_ij.x, index_ij.y - 1), velocityX);
                    float vx2 = checkInside(Index2(index_ij.x - 1, index_ij.y + 1), velocityX);
                    float vx3 = checkInside(Index2(index_ij.x + 1, index_ij.y), velocityX);
                    float vx = bilerp(vx0, vx1, vx2, vx3, 0.5f, 0.5f);

                    Vector2 v(vx, vy);
                    Vector2 prevPoint = point - dt * v;

                    if (prevPoint.x < box.minPosition.x)
                        prevPoint.x = box.minPosition.x;
                    else if (prevPoint.x > box.maxPosition.x)
                        prevPoint.x = box.maxPosition.x;

                    if (prevPoint.y < box.minPosition.y)
                        prevPoint.y = box.minPosition.y;
                    else if (prevPoint.y > box.maxPosition.y)
                        prevPoint.y = box.maxPosition.y;

                    Vector2 indexF = grid.getFaceIndex(prevPoint, 1);
                    Vector2 min(floor(indexF.x), floor(indexF.y));
                    Vector2 max(ceil(indexF.x), ceil(indexF.y));

                    int minX = clamp(min.x, 0, inkRGB.getSize().x - 1);
                    int maxX = clamp(max.x, 0, inkRGB.getSize().x - 1);
                    int minY = clamp(min.y, 0, inkRGB.getSize().y - 1);
                    int maxY = clamp(max.y, 0, inkRGB.getSize().y - 1);

                    Index2 p1(minX, minY);
                    Index2 p2(maxX, minY);
                    Index2 p3(maxX, maxY);
                    Index2 p4(minX, maxY);

                    Vector2 diff = indexF - min;
                    float alpha = diff.x;
                    float betha = diff.y;
                    float vInterpolate = bilerp(velocityY[p1], velocityY[p2], velocityY[p4], velocityY[p3], alpha, betha);

                    VY_result[index_ij] = vInterpolate; 
                }
            }
            velocityX = VX_result;
            velocityY = VY_result;
        }
    }
}

void Fluid2::fluidEmission()
{
    if (Scene::testcase >= Scene::SMOKE) {
        // Emitters contribution HERE
    }
}

void Fluid2::fluidVolumeForces(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) {
        // Gravity term HERE
    }
}

void Fluid2::fluidViscosity(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) {
        // Viscosity term HERE
    }
}

void Fluid2::fluidPressureProjection(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) {
        // Incompressibility / Pressure term HERE
    }
}
}  // namespace asa
