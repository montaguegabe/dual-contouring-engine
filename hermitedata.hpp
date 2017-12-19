//
//  hermitedata.hpp
//  cs175-final
//
//  Created by Gabe Montague on 12/16/17.
//  Copyright © 2017 cs175. All rights reserved.
//

#ifndef hermitedata_hpp
#define hermitedata_hpp

#include "cvec.h"
#include "geometry.h"
#include "dcUtils.h"
#include <set>
#include "implicit.hpp"
#include <memory>
#include <limits>

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
  
  static void pushQuad(std::vector<VertexPNX>& target, const Cvec3f p1, const Cvec3f p2, const Cvec3f p3, const Cvec3f p4, Cvec3f n1, Cvec3f n2, Cvec3f n3, Cvec3f n4, const double scale, const bool flat) {
    
    double diag1Sqr = norm2(p2 - p4);
    double diag2Sqr = norm2(p1 - p3);
    
    VertexPNX v1(p1 * scale, n1, Cvec2f(1, 1));
    VertexPNX v2(p2 * scale, n2, Cvec2f(1, 0));
    VertexPNX v3(p3 * scale, n3, Cvec2f(0, 0));
    VertexPNX v4(p4 * scale, n4, Cvec2f(0, 1));
    
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
    const float thresh = 0.1;
    if (dot(t1v1->n, fNAvg) < thresh || hardEdges) {
      t1v1->n = fN1;
    }
    if (dot(t1v2->n, fNAvg) < thresh || hardEdges) {
      t1v2->n = fN1;
    }
    if (dot(t1v3->n, fNAvg) < thresh || hardEdges) {
      t1v3->n = fN1;
    }
    target.push_back(*t1v1);
    target.push_back(*t1v2);
    target.push_back(*t1v3);
    
    if (dot(t2v1->n, fNAvg) < thresh || hardEdges) {
      t2v1->n = fN2;
    }
    if (dot(t2v2->n, fNAvg) < thresh || hardEdges) {
      t2v2->n = fN2;
    }
    if (dot(t2v3->n, fNAvg) < thresh || hardEdges) {
      t2v3->n = fN2;
    }
    target.push_back(*t2v1);
    target.push_back(*t2v2);
    target.push_back(*t2v3);
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
    
    // Dual contouring randomized
    void triangulateToVector(std::vector<VertexPNX> & target, const double scale, const int randIterations, const bool sharp) {
      
      const int numTrials = randIterations;
      
      for (int x = 0; x < xSize - 1; x++) {
        for (int y = 0; y < ySize - 1; y++) {
          for (int z = 0; z < zSize - 1; z++) {
            
            const Voxel vox = m_edges[x][y][z];
              
            double minError = 1000000;
            Cvec3 minErrorPoint;
            bool setError = false;
            
            
            
            for (int i = 0; i < numTrials; i++) {
              
              const Cvec3 noise(randZeroToOne(), randZeroToOne(), randZeroToOne());
              
              // Calculate sum of squares error
              double error = 0;
              Cvec3 position, diff;
              double dp;
              double dist;
              
              dist = vox.x.dist;
              if (dist >= 0) {
                position = Cvec3(dist, 0, 0);
                diff = noise - position;
                dp = dot(diff, vox.x.normal);
                error += dp * dp;
                setError = true;
              }
              
              dist = vox.y.dist;
              if (dist >= 0) {
                position = Cvec3(0, dist, 0);
                diff = noise - position;
                dp = dot(diff, vox.y.normal);
                error += dp * dp;
                setError = true;
              }
              
              dist = vox.z.dist;
              if (dist >= 0) {
                position = Cvec3(0, 0, dist);
                diff = noise - position;
                dp = dot(diff, vox.z.normal);
                error += dp * dp;
                setError = true;
              }
              
              dist = m_edges[x + 1][y][z].y.dist;
              if (dist >= 0) {
                position = Cvec3(1, dist, 0);
                diff = noise - position;
                dp = dot(diff, m_edges[x + 1][y][z].y.normal);
                error += dp * dp;
                setError = true;
              }
              
              dist = m_edges[x + 1][y][z].z.dist;
              if (dist >= 0) {
                position = Cvec3(1, 0, dist);
                diff = noise - position;
                dp = dot(diff, m_edges[x + 1][y][z].z.normal);
                error += dp * dp;
                setError = true;
              }
              
              dist = m_edges[x][y][z + 1].x.dist;
              if (dist >= 0) {
                position = Cvec3(dist, 0, 1);
                diff = noise - position;
                dp = dot(diff, m_edges[x][y][z + 1].x.normal);
                error += dp * dp;
                setError = true;
              }
              
              dist = m_edges[x][y][z + 1].y.dist;
              if (dist >= 0) {
                position = Cvec3(0, dist, 1);
                diff = noise - position;
                dp = dot(diff, m_edges[x][y][z + 1].y.normal);
                error += dp * dp;
                setError = true;
              }
              
              dist = m_edges[x][y + 1][z].x.dist;
              if (dist >= 0) {
                position = Cvec3(dist, 1, 0);
                diff = noise - position;
                dp = dot(diff, m_edges[x][y + 1][z].x.normal);
                error += dp * dp;
                setError = true;
              }
              
              dist = m_edges[x][y + 1][z].z.dist;
              if (dist >= 0) {
                position = Cvec3(0, 1, dist);
                diff = noise - position;
                dp = dot(diff, m_edges[x][y + 1][z].z.normal);
                error += dp * dp;
                setError = true;
              }
              
              dist = m_edges[x + 1][y][z + 1].y.dist;
              if (dist >= 0) {
                position = Cvec3(1, dist, 1);
                diff = noise - position;
                dp = dot(diff, m_edges[x + 1][y][z + 1].y.normal);
                error += dp * dp;
                setError = true;
              }
              
              dist = m_edges[x + 1][y + 1][z].z.dist;
              if (dist >= 0) {
                position = Cvec3(1, 1, dist);
                diff = noise - position;
                dp = dot(diff, m_edges[x + 1][y + 1][z].z.normal);
                error += dp * dp;
                setError = true;
              }
              
              dist = m_edges[x][y + 1][z + 1].x.dist;
              if (dist >= 0) {
                position = Cvec3(dist, 1, 1);
                diff = noise - position;
                dp = dot(diff, m_edges[x][y + 1][z + 1].x.normal);
                error += dp * dp;
                setError = true;
              }
              
              if (!setError) {
                break;
              }
              
              if (error < minError) {
                minError = error;
                minErrorPoint = noise;
              }
            }
            
            if (setError) {
              // Do something with minimum error point
              Cvec3f minErrorPointF = Cvec3f(minErrorPoint[0], minErrorPoint[1], minErrorPoint[2]);
              m_edges[x][y][z].optPoint = minErrorPointF;
              const Cvec3 norm = m_model.getNorm(Cvec3(x, y, z) + minErrorPoint);
              m_edges[x][y][z].optNormal = Cvec3f(norm[0], norm[1], norm[2]);
              
              /* DEBUG DRAW OF POINTS
              const Cvec3f gVert(x + minErrorPoint[0],
                                 y + minErrorPoint[1],
                                 z + minErrorPoint[2]);
              const Cvec3f p1 = gVert;
              const Cvec3f p2 = gVert + Cvec3f(0.1, 0, 0);
              const Cvec3f p3 = gVert + Cvec3f(0, 0.1, 0.1);
              const VertexPNX v1(p1 * scale, Cvec3f(0, 0, 1), Cvec2f(0, 0));
              const VertexPNX v2(p2 * scale, Cvec3f(0, 0, 1), Cvec2f(1, 0));
              const VertexPNX v3(p3 * scale, Cvec3f(0, 0, 1), Cvec2f(0, 1));
              target.push_back(v1);
              target.push_back(v2);
              target.push_back(v3);
              */
            }
            
          }
        }
      }
      
      // Now connect the found points, cycling through each voxel edge
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
                  const VertexPNX v1(p1 * scale, Cvec3f(0, 0, 1), Cvec2f(0, 0));
                  const VertexPNX v2(p2 * scale, Cvec3f(0, 0, 1), Cvec2f(1, 0));
                  const VertexPNX v3(p3 * scale, Cvec3f(0, 0, 1), Cvec2f(0, 1));
                  target.push_back(v1);
                  target.push_back(v2);
                  target.push_back(v3);
                }
                
                else {
                  const Cvec3f p1(x + xPlane.dist, y, z);
                  const Cvec3f p2(x + xPlane.dist + xPlane.normal[0], y + xPlane.normal[1], z + xPlane.normal[2]);
                  const Cvec3f p3(x + xPlane.dist, y + 0.2f, z);
                  const VertexPNX v1(p1 * scale, Cvec3f(0, 0, 1), Cvec2f(0, 0));
                  const VertexPNX v2(p2 * scale, Cvec3f(0, 0, 1), Cvec2f(1, 0));
                  const VertexPNX v3(p3 * scale, Cvec3f(0, 0, 1), Cvec2f(0, 1));
                  target.push_back(v1);
                  target.push_back(v2);
                  target.push_back(v3);
                }
              }
              
              // Draw y
              const Plane yPlane = edges.y;
              if (yPlane.dist > 0) {
                if (!normals) {
                  const Cvec3f p1(x, y, z);
                  const Cvec3f p2(x + 0.2f, y, z);
                  const Cvec3f p3(x, y + yPlane.dist, z);
                  const VertexPNX v1(p1 * scale, Cvec3f(0, 0, 1), Cvec2f(0, 0));
                  const VertexPNX v2(p2 * scale, Cvec3f(0, 0, 1), Cvec2f(1, 0));
                  const VertexPNX v3(p3 * scale, Cvec3f(0, 0, 1), Cvec2f(0, 1));
                  target.push_back(v1);
                  target.push_back(v2);
                  target.push_back(v3);
                }
                else {
                  const Cvec3f p1(x, y + yPlane.dist, z);
                  const Cvec3f p2(x + 0.2f, y + yPlane.dist, z);
                  const Cvec3f p3(x + yPlane.normal[0], y + yPlane.dist + yPlane.normal[1], z + yPlane.normal[2]);
                  const VertexPNX v1(p1 * scale, Cvec3f(0, 0, 1), Cvec2f(0, 0));
                  const VertexPNX v2(p2 * scale, Cvec3f(0, 0, 1), Cvec2f(1, 0));
                  const VertexPNX v3(p3 * scale, Cvec3f(0, 0, 1), Cvec2f(0, 1));
                  target.push_back(v1);
                  target.push_back(v2);
                  target.push_back(v3);
                }
              }
              
              // Draw z
              const Plane zPlane = edges.z;
              if (zPlane.dist > 0) {
                if (!normals) {
                  const Cvec3f p1(x, y, z);
                  const Cvec3f p2(x, y + 0.2, z);
                  const Cvec3f p3(x, y, z + zPlane.dist);
                  const VertexPNX v1(p1 * scale, Cvec3f(1, 0, 0), Cvec2f(0, 0));
                  const VertexPNX v2(p2 * scale, Cvec3f(1, 0, 0), Cvec2f(1, 0));
                  const VertexPNX v3(p3 * scale, Cvec3f(1, 0, 0), Cvec2f(0, 1));
                  target.push_back(v1);
                  target.push_back(v2);
                  target.push_back(v3);
                } else {
                  const Cvec3f p1(x, y, z + zPlane.dist);
                  const Cvec3f p2(x, y + 0.2f, z + zPlane.dist);
                  const Cvec3f p3(x + zPlane.normal[0], y + zPlane.normal[1], z + zPlane.dist + zPlane.normal[2]);
                  const VertexPNX v1(p1 * scale, Cvec3f(0, 0, 1), Cvec2f(0, 0));
                  const VertexPNX v2(p2 * scale, Cvec3f(0, 0, 1), Cvec2f(1, 0));
                  const VertexPNX v3(p3 * scale, Cvec3f(0, 0, 1), Cvec2f(0, 1));
                  target.push_back(v1);
                  target.push_back(v2);
                  target.push_back(v3);
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
