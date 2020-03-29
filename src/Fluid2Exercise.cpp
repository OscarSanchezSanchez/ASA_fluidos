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
        for (int i = 0; i < grid.getSize().x; i++)
        {
            for (int j = 0; j < grid.getSize().y; j++)
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
                Vector3 inkInterpolate = bilerp(ink[p1], ink[p2], ink[p4], ink[p3], alpha, betha);

                //aqui se deberia comprobar que no se salga
                inkRGB[index_ij] = inkInterpolate;
            }
        }
    }

    {
        // Velocity acvection term HERE
        if (Scene::testcase >= Scene::SMOKE) 
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

                    Vector2 indexPrev = grid.getFaceIndex(prevPoint, 0);
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

                    Vector2 indexPrev = grid.getFaceIndex(prevPoint, 1);
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
        // Emitter Yellow
        AABox2 emitter1(-0.1f, -1.9f, 0.1f, -1.75f);
        Vector3 ink1(1.0f, 1.0f, 0.0f);
        Vector2 vel1(0.0f, 8.0f);

        Vector2 emitterIndexMin1 = grid.getFaceIndex(emitter1.minPosition, 0);
        Vector2 emitterIndexMax1 = grid.getFaceIndex(emitter1.maxPosition, 0);
        Vector2 min1(floor(emitterIndexMin1.x), floor(emitterIndexMin1.y));
        Vector2 max1(ceil(emitterIndexMax1.x), ceil(emitterIndexMax1.y));

        for (int i = min1.x; i < max1.x; i++) {
            for (int j = min1.y; j < max1.y; j++) {
                Index2 ind(i, j);
                inkRGB[ind] = ink1;
                velocityX[ind] = vel1.x;
                velocityY[ind] = vel1.y;
            }
        }

        // Emitter 2
        AABox2 emitter2(-0.2f, -1.9f, -0.11f, -1.75f);
        Vector3 ink2(1.0f, 0.0f, 0.0f);
        Vector2 vel2(0.0f, 8.0f);
 
        Vector2 emitterIndexMin2 = grid.getFaceIndex(emitter2.minPosition, 0);
        Vector2 emitterIndexMax2 = grid.getFaceIndex(emitter2.maxPosition, 0);
        Vector2 min2(floor(emitterIndexMin2.x), floor(emitterIndexMin2.y));
        Vector2 max2(ceil(emitterIndexMax2.x), ceil(emitterIndexMax2.y));

        for (int i = min2.x; i < max2.x; i++) {
            for (int j = min2.y; j < max2.y; j++) {
                Index2 ind(i, j);
                inkRGB[ind] = ink2;
                velocityX[ind] = vel2.x;
                velocityY[ind] = vel2.y;
            }
        }
        // Emitter 3
        AABox2 emitter3(0.11f, -1.9f, 0.2f, -1.75f);
        Vector3 ink3(1.0f, 0.0f, 0.0f);
        Vector2 vel3(0.0f, 8.0f);

        Vector2 emitterIndexMin3 = grid.getFaceIndex(emitter3.minPosition, 0);
        Vector2 emitterIndexMax3 = grid.getFaceIndex(emitter3.maxPosition, 0);
        Vector2 min3(floor(emitterIndexMin3.x), floor(emitterIndexMin3.y));
        Vector2 max3(ceil(emitterIndexMax3.x), ceil(emitterIndexMax3.y));

        for (int i = min3.x; i < max3.x; i++) {
            for (int j = min3.y; j < max3.y; j++) {
                Index2 ind(i, j);
                inkRGB[ind] = ink3;
                velocityX[ind] = vel3.x;
                velocityY[ind] = vel3.y;
            }
        }

    }
}

void Fluid2::fluidVolumeForces(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) {
        // Gravity term HERE
        int size = velocityY.getSize().x * velocityY.getSize().y;
        for (int i = 0; i < size; i++)
        {
            velocityY[i] += dt * Scene::kGravity;
        }
    }
}

void Fluid2::fluidViscosity(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) {
        // Viscosity term HERE
        Vector2 delta = grid.getCellDx();
        float density = Scene::kDensity;
		float viscosity = Scene::kViscosity;

        Index2 faceSizeX = grid.getSizeFacesX();
        Index2 faceSizeY = grid.getSizeFacesY();

        Array2<float> VX_result(velocityX.getSize());
        Array2<float> VY_result(velocityY.getSize());

        for (int i = 0; i < faceSizeX.x; i++)
        {
            for (int j = 0; j < faceSizeX.y; j++)
            {
                Index2 index_ij(i, j);

                float velIJ = checkInside(index_ij, velocityX);

                float velRight = checkInside(Index2(index_ij.x + 1, index_ij.y), velocityX);
                float velTop = checkInside(Index2(index_ij.x, index_ij.y + 1), velocityX);
                float velLeft = checkInside(Index2(index_ij.x - 1, index_ij.y), velocityX);
                float velBottom = checkInside(Index2(index_ij.x, index_ij.y - 1), velocityX);

                float result = velIJ + (dt / density) * viscosity *
                    ((velRight - 2*velIJ  + velLeft)/delta.x + (velTop - 2*velIJ + velBottom)/delta.y);

                VX_result[index_ij] = result; 
            }
        }

        for (int i = 0; i < faceSizeY.x; i++)
        {
            for (int j = 0; j < faceSizeY.y; j++)
            {
                Index2 index_ij(i, j);

                float velIJ = checkInside(index_ij, velocityY);

                float velRight = checkInside(Index2(index_ij.x + 1, index_ij.y), velocityY);
                float velTop = checkInside(Index2(index_ij.x, index_ij.y + 1), velocityY);
                float velLeft = checkInside(Index2(index_ij.x - 1, index_ij.y), velocityY);
                float velBottom = checkInside(Index2(index_ij.x, index_ij.y - 1), velocityY);

                float result = velIJ + (dt / density) * viscosity *
                    ((velRight - 2*velIJ  + velLeft)/delta.x + (velTop - 2*velIJ + velBottom)/delta.y);

                VY_result[index_ij] = result; 
            }
        }
        velocityX = VX_result;
        velocityY = VY_result;        
    }
}

void Fluid2::fluidPressureProjection(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) {
        // Incompressibility / Pressure term HERE
        Index2 gridSize = grid.getSize();

        float density = Scene::kDensity;

        PCGSolver<float> solver;
        SparseMatrix<float> A(gridSize.x * gridSize.y, 5);
        std::vector<float> b(gridSize.x * gridSize.y);
        std::vector<float> X(gridSize.x * gridSize.y);

        Vector2 delta = grid.getCellDx();

        for (int i = 0; i < gridSize.x; i++)
        {
            for (int j = 0; j < gridSize.y; j++)
            {
                //diagonal de la matriz dispersa
                Index2 index_ij(i, j);
                int p_index = pressure.getLinearIndex(i, j);
                if (i == 0 || i == gridSize.x - 1)
                {
                    A.set_element(p_index, p_index, 1.0f / (delta.x * delta.x));
                }
                else
                {
                    A.set_element(p_index, p_index, 2.0f / (delta.x * delta.x));
                }
                //no colisiona con el techo
                if (j == 0)
                {
                    A.add_to_element(p_index, p_index, 1.0f / (delta.y * delta.y));
                }
                else
                {
                    A.add_to_element(p_index, p_index, 2.0f / (delta.y * delta.y));
                }

                //terminos de la matriz
                Index2 p_left(i - 1, j);
				if (p_left.x < gridSize.x)
					A.set_element(p_index, pressure.getLinearIndex(p_left.x, p_left.y), (-1.0f / (delta.x * delta.x)));

				Index2 p_right(i + 1, j);
				if (p_right.x < gridSize.x)
					A.set_element(p_index, pressure.getLinearIndex(p_right.x, p_right.y), (-1.0f / (delta.x * delta.x)));

				Index2 p_bottom(i, j - 1);
				if (p_bottom.y < gridSize.y)
					A.set_element(p_index, pressure.getLinearIndex(p_bottom.x, p_bottom.y), (-1.0f / (delta.y * delta.y)));

				Index2 p_top(i, j + 1);
				if (p_top.y < gridSize.y)
					A.set_element(p_index, pressure.getLinearIndex(p_top.x, p_top.y), (-1.0f / (delta.y * delta.y)));
                
                //rellenamos el vector b
                float vx = checkInside(index_ij, velocityX);
				float vy = checkInside(index_ij, velocityY);
				float vRight = checkInside(p_right, velocityX);
				float vTop = checkInside(p_top, velocityY);

				float deltaVx = (vRight - vx) / delta.x;
				float deltaVy = (vTop - vy) / delta.y;
				float bValue = (-density / dt) * (deltaVx + deltaVy);

				b[p_index] = bValue;
            }
        }
        //resolvemos el sistema
        float residual;
        int its;
        solver.solve(A, b, X, residual, its);

        for (int i = 0; i < gridSize.x; i++)
		{
			for (int j = 0; j < gridSize.y; j++)
			{
				int pres_ind = pressure.getLinearIndex(i, j);
				pressure[Index2(i, j)] = X[pres_ind];
			}
		}

        Array2<float> velX(velocityX.getSize());
		Array2<float> velY(velocityY.getSize());
		for (int i = 0; i < gridSize.x; i++)
		{
			for (int j = 0; j < gridSize.y; j++)
			{
				Index2 index_ij(i, j);

				float ro = checkInside(index_ij, pressure);

				Index2 horizontal = Index2(index_ij.x + 1, index_ij.y);
				float horizontalV = checkInside(horizontal, velocityX);
				float horizontalP = checkInside(horizontal, pressure);
				float resultX = horizontalV - (dt / density) * (horizontalP - ro) / delta.x;

				Index2 szX = velX.getSize();
				velX[horizontal] = resultX;

				Index2 vertical = Index2(index_ij.x, index_ij.y+1);
				float verticalV = checkInside(vertical, velocityY);
				float verticalP = checkInside(vertical, pressure);
				float resultY = verticalV - (dt / density) * (verticalP - ro) / delta.y;

				Index2 szY = velY.getSize();
				velY[vertical] = resultY;
			}
		}
		velocityX = velX;
		velocityY = velY;
    }
}
}  // namespace asa
