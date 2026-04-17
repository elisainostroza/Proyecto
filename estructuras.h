#ifndef ESTRUCTURAS_H
#define ESTRUCTURAS_H


// Estructuras para representar vectores, estado y derivadas
struct Vec {
    double x, y, z;
};

struct Estado {
    Vec r;
    Vec p;
};

struct Derivada {
    Vec drdt;
    Vec dpdt;
};


// Operadores para facilitar operaciones con vectores
inline Vec operator+(const Vec& a, const Vec& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

inline Vec operator-(const Vec& a, const Vec& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

inline Vec operator*(const Vec& v, double s) {
    return {v.x * s, v.y * s, v.z * s};
}

inline Vec operator*(double s, const Vec& v) {
    return v * s;
}

inline Vec cruz(const Vec& a, const Vec& b) {
    return {
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

inline double norma2(const Vec& v) {
    return v.x*v.x + v.y*v.y + v.z*v.z;
}

#endif