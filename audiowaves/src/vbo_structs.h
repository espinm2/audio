#ifndef __VBO_STRUCTS_H__
#define __VBO_STRUCTS_H__

// ======================================================================
// helper structures for VBOs, for rendering (note, the data stored in
// each of these is application specific, adjust as needed!)

struct VBOPosNormalColor {
  VBOPosNormalColor() {}
  VBOPosNormalColor(const glm::vec3 &p) {
    x = p.x; y = p.y; z = p.z;
    nx = 1; ny = 0; nz = 0;
    r = 0; g = 0; b = 0;
  }
  VBOPosNormalColor(const glm::vec3 &p, const glm::vec3 &n, const glm::vec3 &c) {
    x = p.x; y = p.y; z = p.z;
    nx = n.x; ny = n.y; nz = n.z;
    r = c.x; g = c.y; b = c.z;
  }
  float x, y, z;    // position
  float nx, ny, nz; // normal
  float r, g, b;    // color
};

// ======================================================================

struct VBOIndexedTri {
  VBOIndexedTri() {}
  VBOIndexedTri(unsigned int a, unsigned int b, unsigned int c) {
    verts[0] = a;
    verts[1] = b;
    verts[2] = c;
  }
  unsigned int verts[3];
};


// ======================================================================

#endif
