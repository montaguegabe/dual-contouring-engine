//
//  hermitedata.hpp
//  cs175-final
//
//  Created by Gabe Montague on 12/16/17.
//  Copyright © 2017 cs175. All rights reserved.
//

#ifndef hermitedata_hpp
#define hermitedata_hpp

// #define DEBUG_POINTS 1

#include "cvec.h"
#include "geometry.h"
#include "dcUtils.h"
#include <set>
#include "implicit.hpp"
#include <memory>
#include <limits>
#include <algorithm>
#include "Eigen/Dense"

#define BINARY_SEARCH_ITERS 80

namespace DC {
  
  typedef ToggleImplicit ModelType;
  
  // Called "Plane" but represents an edge with Hermite data
  struct Plane {
    Cvec3 normal;
    double dist = -1;
  };
  
  // X, Y, and Z are the 3 edges shooting off in the positive direction
  struct Voxel {
    Voxel() {
    }
    Plane x;
    Plane y;
    Plane z;
    Cvec3f optPoint;
    Cvec3f optNormal;
  };
  
  // Returns relative point of intersection with contour. Binary search
  static inline double edgeContourIntersectionRel(const int x, const int y, const int z, const int direction, const bool isIn /* whether x, y, z is inside */, const ModelType & model) {
    
    const int iterations = BINARY_SEARCH_ITERS;
    
    double upperBound = 1.0;
    double lowerBound = 0.0;
    
    for (int trial = 0; trial < iterations; trial++) {
      const double middle = (lowerBound + upperBound) / 2.0;
      const Cvec3 checkPoint = Cvec3(x, y, z) + Cvec3::direction(direction) * middle;
      bool middleIsIn = model.getValue(checkPoint) < 0.0;
      if (middleIsIn != isIn) {
        upperBound = middle;
      } else {
        lowerBound = middle;
      }
    }
    
    return (upperBound + lowerBound) / 2.0;
  }
  
  static inline void addToLeastSquaresProblem(Eigen::Matrix<double, Eigen::Dynamic, 3> & normals, Eigen::Matrix<double, Eigen::Dynamic, 3> & positions, Eigen::Vector3d & massPoint, const Cvec3 & position, const Cvec3 & normal, bool & marked) {
    
    using namespace Eigen;
    const int rows = normals.rows();
    normals.conservativeResize(rows + 1, NoChange);
    positions.conservativeResize(rows + 1, NoChange);
    normals(rows, 0) = normal[0];
    normals(rows, 1) = normal[1];
    normals(rows, 2) = normal[2];
    massPoint += Vector3d({position[0], position[1], position[2]});
    positions(rows, 0) = position[0];
    positions(rows, 1) = position[1];
    positions(rows, 2) = position[2];
    //positions(rows) = dot(normal, position);
  }
  
  template <int xSize, int ySize, int zSize>
  class HermiteData {
  public:
    HermiteData(ModelType model) : m_model(model) {
      
      for (int x = 0; x < xSize - 1; x++) {
        for (int y = 0; y < ySize - 1; y++) {
          for (int z = 0; z < zSize - 1; z++) {
            
            const Voxel zero;
            
            const Cvec3 pos(x, y, z);
            bool isIn = m_model.getValue(pos) < 0.0;
            
            // Determine x extent
            if ((m_model.getValue(Cvec3(x + 1, y, z)) < 0.0) != isIn) {
              
              const double distance = edgeContourIntersectionRel(x, y, z, 0, isIn, m_model);
              Plane & targetEdge = m_edges[x][y][z].x;
              targetEdge.dist = distance;
              targetEdge.normal = m_model.getNorm(Cvec3(x + distance, y, z));
            }
            
            // Determine y extent
            if ((m_model.getValue(Cvec3(x, y + 1, z)) < 0.0) != isIn) {
              const double distance = edgeContourIntersectionRel(x, y, z, 1, isIn, m_model);
              Plane & targetEdge = m_edges[x][y][z].y;
              targetEdge.dist = distance;
              targetEdge.normal = m_model.getNorm(Cvec3(x, y + distance, z));
            }
            
            // Determine z extent
            if ((m_model.getValue(Cvec3(x, y, z + 1)) < 0.0) != isIn) {
              const double distance = edgeContourIntersectionRel(x, y, z, 2, isIn, m_model);
              Plane & targetEdge = m_edges[x][y][z].z;
              targetEdge.dist = distance;
              targetEdge.normal = m_model.getNorm(Cvec3(x, y, z + distance));
            }
          }
        }
      }
    }
    
    HermiteData() {
      
      for (int x = 0; x < xSize - 1; x++) {
        for (int y = 0; y < ySize - 1; y++) {
          for (int z = 0; z < zSize - 1; z++) {
            
            const Voxel zero;
            m_edges[x][y][z] = zero;
          }
        }
      }
    }
    
    // Members
    ModelType m_model;
    Voxel m_edges[xSize][ySize][zSize];
    
    static void writeVertex(std::vector<VertexPNX> & target, double scale, VertexPNX v) {
      v.p = (v.p - Cvec3f(0.5 * xSize, 0.5 * ySize, 0.5 * zSize)) * scale;
      target.push_back(v);
    }
    
    static void pushQuad(std::vector<VertexPNX>& target, const Cvec3f p1, const Cvec3f p2, const Cvec3f p3, const Cvec3f p4, Cvec3f n1, Cvec3f n2, Cvec3f n3, Cvec3f n4, const double scale, const bool flat) {
      
      double diag1Sqr = norm2(p2 - p4);
      double diag2Sqr = norm2(p1 - p3);
      
      VertexPNX v1(p1, n1, Cvec2f(1, 1));
      VertexPNX v2(p2, n2, Cvec2f(1, 0));
      VertexPNX v3(p3, n3, Cvec2f(0, 0));
      VertexPNX v4(p4, n4, Cvec2f(0, 1));
      
      VertexPNX *t1v1, *t1v2, *t1v3, *t2v1, *t2v2, *t2v3;
      const Cvec3f *t1p1, *t1p2, *t1p3, *t2p1, *t2p2, *t2p3;
      if (diag1Sqr > diag2Sqr) {
        t1v1 = &v3;
        t1p1 = &p3;
        t1v2 = &v2;
        t1p2 = &p2;
        t1v3 = &v1;
        t1p3 = &p1;
        
        t2v1 = &v3;
        t2p1 = &p3;
        t2v2 = &v1;
        t2p2 = &p1;
        t2v3 = &v4;
        t2p3 = &p4;
        
      } else {
        // 3 -> 4 -> 1 -> 2 ->
        t1v1 = &v4;
        t1p1 = &p4;
        t1v2 = &v3;
        t1p2 = &p3;
        t1v3 = &v2;
        t1p3 = &p2;
        
        t2v1 = &v4;
        t2p1 = &p4;
        t2v2 = &v2;
        t2p2 = &p2;
        t2v3 = &v1;
        t2p3 = &p1;
      }
      
      const Cvec3f nAvg = n1 + n2 + n3 + n4;
      Cvec3f fN1 = normalize(cross(*t1p3 - *t1p1, *t1p2 - *t1p1));
      if (dot(fN1, nAvg) < 0.0f) {
        fN1 = fN1 * -1.0f;
      }
      Cvec3f fN2 = normalize(cross(*t2p3 - *t2p1, *t2p2 - *t2p1));
      if (dot(fN2, nAvg) < 0.0f) {
        fN2 = fN2 * -1.0f;
      }
      Cvec3f fNAvg = normalize(fN1 + fN2);
      
      bool hardEdges = flat;
      const float thresh = -100000;
      if (dot(t1v1->n, fNAvg) < thresh || hardEdges) {
        t1v1->n = fN1;
      }
      if (dot(t1v2->n, fNAvg) < thresh || hardEdges) {
        t1v2->n = fN1;
      }
      if (dot(t1v3->n, fNAvg) < thresh || hardEdges) {
        t1v3->n = fN1;
      }
      writeVertex(target, scale, *t1v1);
      writeVertex(target, scale, *t1v2);
      writeVertex(target, scale, *t1v3);
      
      if (dot(t2v1->n, fNAvg) < thresh || hardEdges) {
        t2v1->n = fN2;
      }
      if (dot(t2v2->n, fNAvg) < thresh || hardEdges) {
        t2v2->n = fN2;
      }
      if (dot(t2v3->n, fNAvg) < thresh || hardEdges) {
        t2v3->n = fN2;
      }
      writeVertex(target, scale, *t2v1);
      writeVertex(target, scale, *t2v2);
      writeVertex(target, scale, *t2v3);
    }
    
    // Dual contouring
    void triangulateToVector(std::vector<VertexPNX> & target, const double scale, const double tolerance, const bool sharp) {
      
      using std::vector;
      using namespace Eigen;
      
      for (int x = 0; x < xSize - 1; x++) {
        for (int y = 0; y < ySize - 1; y++) {
          for (int z = 0; z < zSize - 1; z++) {
            
            const Voxel vox = m_edges[x][y][z];
            
            // Calculate sum of squares error
            Matrix<double, Dynamic, 3> normals;
            Matrix<double, Dynamic, 3> positions;
            Vector3d massPoint = Vector3d::Zero();
            
            double dist;
            Cvec3 position;
            Cvec3 normal;
            bool marked = false;
            
            dist = vox.x.dist;
            if (dist >= 0) {
              position = Cvec3(dist, 0, 0);
              normal = vox.x.normal;
              addToLeastSquaresProblem(normals, positions, massPoint, position, normal, marked);
            }
            
            dist = vox.y.dist;
            if (dist >= 0) {
              position = Cvec3(0, dist, 0);
              normal = vox.y.normal;
              addToLeastSquaresProblem(normals, positions, massPoint, position, normal, marked);
            }
            
            dist = vox.z.dist;
            if (dist >= 0) {
              position = Cvec3(0, 0, dist);
              normal = vox.z.normal;
              addToLeastSquaresProblem(normals, positions, massPoint, position, normal, marked);
            }
            
            dist = m_edges[x + 1][y][z].y.dist;
            if (dist >= 0) {
              position = Cvec3(1, dist, 0);
              normal = m_edges[x + 1][y][z].y.normal;
              addToLeastSquaresProblem(normals, positions, massPoint, position, normal, marked);
            }
            
            dist = m_edges[x + 1][y][z].z.dist;
            if (dist >= 0) {
              position = Cvec3(1, 0, dist);
              normal = m_edges[x + 1][y][z].z.normal;
              addToLeastSquaresProblem(normals, positions, massPoint, position, normal, marked);
            }
            
            dist = m_edges[x][y][z + 1].x.dist;
            if (dist >= 0) {
              position = Cvec3(dist, 0, 1);
              normal = m_edges[x][y][z + 1].x.normal;
              addToLeastSquaresProblem(normals, positions, massPoint, position, normal, marked);
            }
            
            dist = m_edges[x][y][z + 1].y.dist;
            if (dist >= 0) {
              position = Cvec3(0, dist, 1);
              normal = m_edges[x][y][z + 1].y.normal;
              addToLeastSquaresProblem(normals, positions, massPoint, position, normal, marked);
            }
            
            dist = m_edges[x][y + 1][z].x.dist;
            if (dist >= 0) {
              position = Cvec3(dist, 1, 0);
              normal = m_edges[x][y + 1][z].x.normal;
              addToLeastSquaresProblem(normals, positions, massPoint, position, normal, marked);
            }
            
            dist = m_edges[x][y + 1][z].z.dist;
            if (dist >= 0) {
              position = Cvec3(0, 1, dist);
              normal = m_edges[x][y + 1][z].z.normal;
              addToLeastSquaresProblem(normals, positions, massPoint, position, normal, marked);
            }
            
            dist = m_edges[x + 1][y][z + 1].y.dist;
            if (dist >= 0) {
              position = Cvec3(1, dist, 1);
              normal = m_edges[x + 1][y][z + 1].y.normal;
              addToLeastSquaresProblem(normals, positions, massPoint, position, normal, marked);
            }
            
            dist = m_edges[x + 1][y + 1][z].z.dist;
            if (dist >= 0) {
              position = Cvec3(1, 1, dist);
              normal = m_edges[x + 1][y + 1][z].z.normal;
              addToLeastSquaresProblem(normals, positions, massPoint, position, normal, marked);
            }
            
            dist = m_edges[x][y + 1][z + 1].x.dist;
            if (dist >= 0) {
              position = Cvec3(dist, 1, 1);
              normal = m_edges[x][y + 1][z + 1].x.normal;
              addToLeastSquaresProblem(normals, positions, massPoint, position, normal, marked);
            }
            
            const int rows = normals.rows();
            if (rows > 0) {
              
              using std::cout;
              using std::endl;
              
              
              // Step 0: Calculate mass point
              massPoint = massPoint / rows;
              
              // Calculate real positions values
              Matrix<double, Dynamic, 1> ideals;
              ideals.conservativeResize(rows, NoChange);
              for (int i = 0; i < rows; i++) {
                ideals(i) = (
                   normals(i, 0) * (positions(i, 0) - massPoint(0)) +
                   normals(i, 1) * (positions(i, 1) - massPoint(1)) +
                   normals(i, 2) * (positions(i, 2) - massPoint(2))
                                 );
              }
              
              // Step 1: Find QR decomposition of merged
              Matrix<double, Dynamic, 4> planeMatrix(normals.rows(), 4);
              planeMatrix << normals, ideals;
              const HouseholderQR<Matrix<double, Dynamic, 4> > decomp = planeMatrix.householderQr();
              Matrix<double, Dynamic, Dynamic> matrixQ = decomp.householderQ();
              Matrix<double, Dynamic, 4> matrixR = (decomp.matrixQR().triangularView<Upper>());
              
              if (rows < 3) {
                matrixR.conservativeResizeLike(MatrixXd::Zero(3, 4));
              }
              
              const Matrix<double, 3, 3> matrixAHat = matrixR.block(0, 0, 3, 3);
              const Matrix<double, 3, 1> vectorBHat = matrixR.block(0, 3, 3, 1);
              
              // Step 2: Compute the SVD decomposition of A hat
              /*const Matrix<double, 3, 3> matrixAT = matrixAHat.transpose();
              const Matrix<double, 3, 3> matrixATA = matrixAT * matrixAHat;
              //const Matrix<double, 3, 3> matrixAAT = matrixAHat * matrixAT;
              
              SelfAdjointEigenSolver<Matrix<double, 3, 3> > eigensolver(matrixATA);
              Vector3d eigenvalues = eigensolver.eigenvalues();
              Matrix<double, 3, 3> matrixDInv;
              Matrix<double, 3, 3> matrixD;
              const double tolerance = 0.000001;
              for (int i = 0; i < 3; i++) {
                const double singularValue = std::sqrt(eigenvalues(i));
                matrixDInv(i, i) = singularValue < tolerance ? 0.0 : 1.0 / singularValue;
                matrixD(i, i) = singularValue;
                if (singularValue < 0.1) {
                  //marked = true;
                }
              }*/
              
              // Rows of U are eigenvectors
              /*const Matrix<double, 3, 3> matrixUT = eigensolver.eigenvectors();
              const Matrix<double, 3, 3> matrixU = matrixUT.transpose();
              const Matrix<double, Dynamic, Dynamic> svdProduct = matrixUT * matrixD * matrixU;
              const Matrix<double, Dynamic, Dynamic> svdProduct2 = matrixU * matrixD * matrixUT;
              
              cout << "---" << endl;
              cout << matrixAHat << endl;
              cout << "---" << endl;
              cout << svdProduct << endl;
              cout << "---" << endl;
              cout << svdProduct2 << endl;
              cout << "---" << endl;
              cout << "---" << endl;
              */
              JacobiSVD<Matrix<double, 3, 3> > svd(matrixAHat, ComputeFullU | ComputeFullV);
              Matrix<double, 3, 3> pseudoInv;
              svd.pinv(pseudoInv, tolerance);
              
              //const Matrix<double, 3, 1> solution =  matrixAHat.jacobiSvd(ComputeFullU | ComputeFullV).solve(vectorBHat);
              const Matrix<double, 3, 1> solution = pseudoInv * vectorBHat;
 
              // Add solution point to structure
              Cvec3 minErrorPoint = Cvec3(solution(0) + massPoint(0), solution(1) + massPoint(1), solution(2) + massPoint(2));
              if (rows < 3 || marked) {
                // Underdetermined
                minErrorPoint[0] = 0;
                minErrorPoint[1] = -10;
                minErrorPoint[2] = 0;
              }
              Cvec3f minErrorPointF = Cvec3f(minErrorPoint[0], minErrorPoint[1], minErrorPoint[2]);
              m_edges[x][y][z].optPoint = minErrorPointF;
              const Cvec3 norm = m_model.getNorm(Cvec3(x, y, z) + minErrorPoint);
              m_edges[x][y][z].optNormal = Cvec3f(norm[0], norm[1], norm[2]);
              
              // DEBUG DRAW OF POINTS
#ifdef DEBUG_POINTS
              const Cvec3f gVert(x + minErrorPoint[0],
                                 y + minErrorPoint[1],
                                 z + minErrorPoint[2]);
              const Cvec3f p1 = gVert;
              const Cvec3f p2 = gVert + Cvec3f(0.1, 0, 0);
              const Cvec3f p3 = gVert + Cvec3f(0, 0.1, 0.1);
              const VertexPNX v1(p1, Cvec3f(0, 0, 1), Cvec2f(0, 0));
              const VertexPNX v2(p2, Cvec3f(0, 0, 1), Cvec2f(1, 0));
              const VertexPNX v3(p3, Cvec3f(0, 0, 1), Cvec2f(0, 1));
              writeVertex(target, scale, v1);
              writeVertex(target, scale, v2);
              writeVertex(target, scale, v3);
#endif
            }
            
          }
        }
      }
      
      // Now connect the found points, cycling through each voxel edge
#ifndef DEBUG_POINTS
      for (int x = 1; x < xSize - 1; x++) {
        for (int y = 1; y < ySize - 1; y++) {
          for (int z = 1; z < zSize - 1; z++) {
            
            if (m_edges[x][y][z].x.dist > 0) {
              const Cvec3f p1 = Cvec3f(x, y, z) + m_edges[x][y][z].optPoint;
              Cvec3f n1 = m_edges[x][y][z].optNormal;
              const Cvec3f p2 = Cvec3f(x, y - 1, z) + m_edges[x][y-1][z].optPoint;
              Cvec3f n2 = m_edges[x][y - 1][z].optNormal;
              const Cvec3f p3 = Cvec3f(x, y - 1, z - 1) + m_edges[x][y-1][z-1].optPoint;
              Cvec3f n3 = m_edges[x][y - 1][z - 1].optNormal;
              const Cvec3f p4 = Cvec3f(x, y, z - 1) + m_edges[x][y][z-1].optPoint;
              Cvec3f n4 = m_edges[x][y][z - 1].optNormal;
              pushQuad(target, p1, p2, p3, p4, n1, n2, n3, n4, scale, sharp);
            }
            
            if (m_edges[x][y][z].y.dist > 0) {
              const Cvec3f p1 = Cvec3f(x, y, z) + m_edges[x][y][z].optPoint;
              Cvec3f n1 = m_edges[x][y][z].optNormal;
              const Cvec3f p2 = Cvec3f(x - 1, y, z) + m_edges[x-1][y][z].optPoint;
              Cvec3f n2 = m_edges[x - 1][y][z].optNormal;
              const Cvec3f p3 = Cvec3f(x - 1, y, z - 1) + m_edges[x-1][y][z-1].optPoint;
              Cvec3f n3 = m_edges[x - 1][y][z - 1].optNormal;
              const Cvec3f p4 = Cvec3f(x, y, z - 1) + m_edges[x][y][z-1].optPoint;
              Cvec3f n4 = m_edges[x][y][z - 1].optNormal;
              pushQuad(target, p1, p2, p3, p4, n1, n2, n3, n4, scale, sharp);
            }
            
            if (m_edges[x][y][z].z.dist > 0) {
              const Cvec3f p1 = Cvec3f(x, y, z) + m_edges[x][y][z].optPoint;
              Cvec3f n1 = m_edges[x][y][z].optNormal;
              const Cvec3f p2 = Cvec3f(x - 1, y, z) + m_edges[x-1][y][z].optPoint;
              Cvec3f n2 = m_edges[x - 1][y][z].optNormal;
              const Cvec3f p3 = Cvec3f(x - 1, y - 1, z) + m_edges[x-1][y-1][z].optPoint;
              Cvec3f n3 = m_edges[x - 1][y - 1][z].optNormal;
              const Cvec3f p4 = Cvec3f(x, y - 1, z) + m_edges[x][y-1][z].optPoint;
              Cvec3f n4 = m_edges[x][y - 1][z].optNormal;
              pushQuad(target, p1, p2, p3, p4, n1, n2, n3, n4, scale, sharp);
            }
          }
        }
      }
#endif
    }
    
    // DEBUG function to display distances and normals
    void triangulateToVectorDebug(std::vector<VertexPNX> & target, const double scale, const bool normals) const {
      
      for (int x = 0; x < xSize; x++) {
        for (int y = 0; y < ySize; y++) {
          for (int z = 0; z < zSize; z++) {
      
            const Voxel edges = m_edges[x][y][z];
            if (true) {
              
              // Draw x
              const Plane xPlane = edges.x;
              if (xPlane.dist > 0) {
                if (!normals) {
                  const Cvec3f p1(x, y, z);
                  const Cvec3f p2(x + xPlane.dist, y, z);
                  const Cvec3f p3(x, y + 0.2f, z);
                  const VertexPNX v1(p1, Cvec3f(0, 0, 1), Cvec2f(0, 0));
                  const VertexPNX v2(p2, Cvec3f(0, 0, 1), Cvec2f(1, 0));
                  const VertexPNX v3(p3, Cvec3f(0, 0, 1), Cvec2f(0, 1));
                  writeVertex(target, scale, v1);
                  writeVertex(target, scale, v2);
                  writeVertex(target, scale, v3);
                }
                
                else {
                  const Cvec3f p1(x + xPlane.dist, y, z);
                  const Cvec3f p2(x + xPlane.dist + xPlane.normal[0], y + xPlane.normal[1], z + xPlane.normal[2]);
                  const Cvec3f p3(x + xPlane.dist, y + 0.2f, z);
                  const VertexPNX v1(p1, Cvec3f(0, 0, 1), Cvec2f(0, 0));
                  const VertexPNX v2(p2, Cvec3f(0, 0, 1), Cvec2f(1, 0));
                  const VertexPNX v3(p3, Cvec3f(0, 0, 1), Cvec2f(0, 1));
                  writeVertex(target, scale, v1);
                  writeVertex(target, scale, v2);
                  writeVertex(target, scale, v3);
                }
              }
              
              // Draw y
              const Plane yPlane = edges.y;
              if (yPlane.dist > 0) {
                if (!normals) {
                  const Cvec3f p1(x, y, z);
                  const Cvec3f p2(x + 0.2f, y, z);
                  const Cvec3f p3(x, y + yPlane.dist, z);
                  const VertexPNX v1(p1, Cvec3f(0, 0, 1), Cvec2f(0, 0));
                  const VertexPNX v2(p2, Cvec3f(0, 0, 1), Cvec2f(1, 0));
                  const VertexPNX v3(p3, Cvec3f(0, 0, 1), Cvec2f(0, 1));
                  writeVertex(target, scale, v1);
                  writeVertex(target, scale, v2);
                  writeVertex(target, scale, v3);
                }
                else {
                  const Cvec3f p1(x, y + yPlane.dist, z);
                  const Cvec3f p2(x + 0.2f, y + yPlane.dist, z);
                  const Cvec3f p3(x + yPlane.normal[0], y + yPlane.dist + yPlane.normal[1], z + yPlane.normal[2]);
                  const VertexPNX v1(p1, Cvec3f(0, 0, 1), Cvec2f(0, 0));
                  const VertexPNX v2(p2, Cvec3f(0, 0, 1), Cvec2f(1, 0));
                  const VertexPNX v3(p3, Cvec3f(0, 0, 1), Cvec2f(0, 1));
                  writeVertex(target, scale, v1);
                  writeVertex(target, scale, v2);
                  writeVertex(target, scale, v3);
                }
              }
              
              // Draw z
              const Plane zPlane = edges.z;
              if (zPlane.dist > 0) {
                if (!normals) {
                  const Cvec3f p1(x, y, z);
                  const Cvec3f p2(x, y + 0.2, z);
                  const Cvec3f p3(x, y, z + zPlane.dist);
                  const VertexPNX v1(p1, Cvec3f(1, 0, 0), Cvec2f(0, 0));
                  const VertexPNX v2(p2, Cvec3f(1, 0, 0), Cvec2f(1, 0));
                  const VertexPNX v3(p3, Cvec3f(1, 0, 0), Cvec2f(0, 1));
                  writeVertex(target, scale, v1);
                  writeVertex(target, scale, v2);
                  writeVertex(target, scale, v3);
                } else {
                  const Cvec3f p1(x, y, z + zPlane.dist);
                  const Cvec3f p2(x, y + 0.2f, z + zPlane.dist);
                  const Cvec3f p3(x + zPlane.normal[0], y + zPlane.normal[1], z + zPlane.dist + zPlane.normal[2]);
                  const VertexPNX v1(p1, Cvec3f(0, 0, 1), Cvec2f(0, 0));
                  const VertexPNX v2(p2, Cvec3f(0, 0, 1), Cvec2f(1, 0));
                  const VertexPNX v3(p3, Cvec3f(0, 0, 1), Cvec2f(0, 1));
                  writeVertex(target, scale, v1);
                  writeVertex(target, scale, v2);
                  writeVertex(target, scale, v3);
                }
                
              }
            }

          }
        }
      }
    }
  };
  
  
}

#endif /* hermitedata_hpp */
