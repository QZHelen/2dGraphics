#include "GBitmap.h"
#include "GCanvas.h"
#include "GColor.h"
#include "GMath.h"
#include "GPixel.h"
#include "GPoint.h"
#include "GPaint.h"
#include "GRect.h"
#include <algorithm>
#include <vector>

#include <iostream>

#define BYTE_MAX_VALUE  255
#define HALF_BYTE_VALUE  127
/** PA 2 */
struct Edge { // edge struct for polygons
    int top_y;
    int bottom_y;
    float curr_x;
    float slope;

    bool operator < (const Edge& e) const
    {
        return (top_y < e.top_y) | 
            ((top_y == e.top_y) & (curr_x < e.curr_x)) | 
            ((top_y == e.top_y) & (curr_x == e.curr_x) & (slope < e.slope));
    }
} ;
/**
 *  Given 2 points, add valid edge(s) made from them. 0 <= #edge <= 3 to the end of 
 *  Edge vector parameter.
 */
static inline void clipPolygonPoints(
    GPoint p0, GPoint p1, float top, float bottom, float left, float right, std::vector<Edge> & edges) {
    // (1) Horizontal case prep: sort the 2 points by y value, small to large.
    if (p0.y() > p1.y()) {
        GPoint temp = p0;
        p0 = p1;
        p1 = temp;
    }
    
    // calculate slope, terminate function if both y's are equal.
    int top_y = GRoundToInt(p0.y());
    int bottom_y = GRoundToInt(p1.y());
    if (top_y == bottom_y) {
        return;
    }
    float slope = (p1.x() - p0.x()) / (p1.y() - p0.y());
    
    // completely above or below the bounds: exit.
    if (((p1.y() <= top) | (p0.y() >= bottom))) {
        return;
    }
    // one point above or below the bounds: trim and calculating new x using 
    // newX = slope * (newY - p0.y) + p0.x
    if (p0.y() < top) {
       p0.set(slope * (top - p0.y()) + p0.x(), top); 
    }
    if (p1.y() > bottom) {
        p1.set(slope * (bottom - p0.y()) + p0.x(), bottom);
    }
    // recalculate top and bottom y for edges
    top_y = GRoundToInt(p0.y());
    bottom_y = GRoundToInt(p1.y());

    // (2) Vertical case prep: sort the 2 points by x value, small to large.
    if (p0.x() > p1.x()) {
        GPoint temp = p0;
        p0 = p1;
        p1 = temp;
    }
    // completely to the left or right of the bounds: project both X'es.
    if (p1.x() <= left) {
        p0.set(left, p0.y());
        p1.set(left, p1.y());
        slope = 0.0f;
        // edges.push_back({top_y, bottom_y, left, 0.0f});
        // return;
    }
    if (p0.x() >= right) {
        p0.set(right, p0.y());
        p1.set(right, p1.y());
        slope = 0.0f;
        // edges.push_back({top_y, bottom_y, right, 0.0f});
        // return;
    }
    // only one to the left or right: chop and project. calculating y intersection using
    // newY = (newX - p0.x) / slope + p0.y
    if (p0.x() < left) {
        float intersectY = (left - p0.x()) / slope + p0.y();
        top_y =  GRoundToInt(std::min(intersectY, p0.y()));
        bottom_y =  GRoundToInt(std::max(intersectY, p0.y()));
        edges.push_back({top_y, bottom_y, left, 0.0f});
        p0.set(left, intersectY);
    }
    if (p1.x() > right) {
        float intersectY = (right - p1.x()) / slope + p1.y();
        top_y =  GRoundToInt(std::min(intersectY, p1.y()));
        bottom_y = GRoundToInt(std::max(intersectY, p1.y()));
        edges.push_back({top_y, bottom_y, left, 0.0f});
        p1.set(right, intersectY);
    }
    // resort p0, p1 by y and add the in-bound edge.
    if (p0.y() > p1.y()) {
        GPoint temp = p0;
        p0 = p1;
        p1 = temp;
    }
    top_y =  GRoundToInt(p0.y());
    bottom_y =  GRoundToInt(p1.y());
    if (top_y != bottom_y) {
        edges.push_back({top_y, bottom_y, p0.x(), slope});
    }
}
/** PA 1-2 */
/**
 *  Given a GColor with ARGB values in range [0, 1], convert it to a GPixel in range [0, 255].
 */
static inline GPixel convertColorToPixel(const GColor& c) {
    // pin to unit.
    GColor safeC = c.pinToUnit();
    float a = safeC.fA;
    float r = safeC.fR;
    float g = safeC.fG;
    float b = safeC.fB;
    
    // premul rgb.
    r *= a;
    g *= a;
    b *= a;
    
    // scale to 255 with rounding.
    a = GRoundToInt(a * BYTE_MAX_VALUE);
    r = GRoundToInt(r * BYTE_MAX_VALUE);
    g = GRoundToInt(g * BYTE_MAX_VALUE);
    b = GRoundToInt(b * BYTE_MAX_VALUE);
    
    // pack into GPixel.
    return GPixel_PackARGB(a, r, g, b);
}

/**
 *  Return the result of multiplying 2 bytes [0, 255] values and return a byte value [0, 255].
 */
static inline int multiplyBytes(int b1, int b2) {
    int r = b1 * b2;
    return (r + HALF_BYTE_VALUE) / BYTE_MAX_VALUE; /* re-normalize and round the result. */
}

/**
 *  Below are 12 helper functions for calculating the result pixel for each of 
 *  the 12 blend modes. 
 */
static GPixel kClear(const GPixel source, const GPixel dest) {
    //!< [0, 0]
    return GPixel_PackARGB(0, 0, 0, 0);
}
static GPixel kSrc(const GPixel source, const GPixel dest) {
    // TODO: should make a copy here?
    //!< [Sa, Sc]
    return source;
}
static GPixel kDst(const GPixel source, const GPixel dest) {
    //!< [Da, Dc]
    return dest;
}
static GPixel kSrcOver(const GPixel source, const GPixel dest) {
    //!< [Sa + Da * (1 - Sa), Sc + Dc * (1 - Sa)]
    int sA = GPixel_GetA(source);
    int sR = GPixel_GetR(source);
    int sG = GPixel_GetG(source);
    int sB = GPixel_GetB(source);
    
    int dA = GPixel_GetA(dest);
    int dR = GPixel_GetR(dest);
    int dG = GPixel_GetG(dest);
    int dB = GPixel_GetB(dest);
    
    int resultA = sA + multiplyBytes(BYTE_MAX_VALUE - sA, dA);
    int resultR = sR + multiplyBytes(BYTE_MAX_VALUE - sA, dR);
    int resultG = sG + multiplyBytes(BYTE_MAX_VALUE - sA, dG);
    int resultB = sB + multiplyBytes(BYTE_MAX_VALUE - sA, dB);   
    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}
static GPixel kDstOver(const GPixel source, const GPixel dest) {
    //!< [Da + Sa * (1 - Da), Dc + Sc * (1 - Da)]
    int sA = GPixel_GetA(source);
    int sR = GPixel_GetR(source);
    int sG = GPixel_GetG(source);
    int sB = GPixel_GetB(source);
    
    int dA = GPixel_GetA(dest);
    int dR = GPixel_GetR(dest);
    int dG = GPixel_GetG(dest);
    int dB = GPixel_GetB(dest);
    
    int resultA = dA + multiplyBytes(BYTE_MAX_VALUE - dA, sA);
    int resultR = dR + multiplyBytes(BYTE_MAX_VALUE - dA, sR);
    int resultG = dG + multiplyBytes(BYTE_MAX_VALUE - dA, sG);
    int resultB = dB + multiplyBytes(BYTE_MAX_VALUE - dA, sB); 
    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}
static GPixel kSrcIn(const GPixel source, const GPixel dest) {
    //!< [Sa * Da, Sc * Da]
    int sA = GPixel_GetA(source);
    int sR = GPixel_GetR(source);
    int sG = GPixel_GetG(source);
    int sB = GPixel_GetB(source);
    
    int dA = GPixel_GetA(dest);

    int resultA = multiplyBytes(sA, dA);
    int resultR = multiplyBytes(sR, dA);
    int resultG = multiplyBytes(sG, dA);
    int resultB = multiplyBytes(sB, dA);
    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}
static GPixel kDstIn(const GPixel source, const GPixel dest) {
    //!< [Da * Sa, Dc * Sa]  
    int sA = GPixel_GetA(source);
    
    int dA = GPixel_GetA(dest);
    int dR = GPixel_GetR(dest);
    int dG = GPixel_GetG(dest);
    int dB = GPixel_GetB(dest);

    int resultA = multiplyBytes(sA, dA);
    int resultR = multiplyBytes(sA, dR);
    int resultG = multiplyBytes(sA, dG);
    int resultB = multiplyBytes(sA, dB);
    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}
static GPixel kSrcOut(const GPixel source, const GPixel dest) {
    //!< [Sa * (1 - Da), Sc * (1 - Da)]
    int sA = GPixel_GetA(source);
    int sR = GPixel_GetR(source);
    int sG = GPixel_GetG(source);
    int sB = GPixel_GetB(source);
    
    int dA = GPixel_GetA(dest);

    int resultA = multiplyBytes(sA, BYTE_MAX_VALUE - dA);
    int resultR = multiplyBytes(sR, BYTE_MAX_VALUE - dA);
    int resultG = multiplyBytes(sG, BYTE_MAX_VALUE - dA);
    int resultB = multiplyBytes(sB, BYTE_MAX_VALUE - dA);
    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}
static GPixel kDstOut(const GPixel source, const GPixel dest) {
    //!< [Da * (1 - Sa), Dc * (1 - Sa)]
    int sA = GPixel_GetA(source);
    
    int dA = GPixel_GetA(dest);
    int dR = GPixel_GetR(dest);
    int dG = GPixel_GetG(dest);
    int dB = GPixel_GetB(dest);

    int resultA = multiplyBytes(dA, BYTE_MAX_VALUE - sA);
    int resultR = multiplyBytes(dR, BYTE_MAX_VALUE - sA);
    int resultG = multiplyBytes(dG, BYTE_MAX_VALUE - sA);
    int resultB = multiplyBytes(dB, BYTE_MAX_VALUE - sA);
    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}
static GPixel kSrcATop(const GPixel source, const GPixel dest) {
    //!< [Da, Sc * Da + Dc * (1 - Sa)]
    int sA = GPixel_GetA(source);
    int sR = GPixel_GetR(source);
    int sG = GPixel_GetG(source);
    int sB = GPixel_GetB(source);
    
    int dA = GPixel_GetA(dest);
    int dR = GPixel_GetR(dest);
    int dG = GPixel_GetG(dest);
    int dB = GPixel_GetB(dest);

    int resultA = dA;
    int resultR = multiplyBytes(sR, dA) + multiplyBytes(dR, BYTE_MAX_VALUE - sA);
    int resultG = multiplyBytes(sG, dA) + multiplyBytes(dG, BYTE_MAX_VALUE - sA);
    int resultB = multiplyBytes(sB, dA) + multiplyBytes(dB, BYTE_MAX_VALUE - sA);
    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}
static GPixel kDstATop(const GPixel source, const GPixel dest) {
    //!< [Sa, Dc * Sa + Sc * (1 - Da)]
    int sA = GPixel_GetA(source);
    int sR = GPixel_GetR(source);
    int sG = GPixel_GetG(source);
    int sB = GPixel_GetB(source);
    
    int dA = GPixel_GetA(dest);
    int dR = GPixel_GetR(dest);
    int dG = GPixel_GetG(dest);
    int dB = GPixel_GetB(dest);

    int resultA = sA;
    int resultR = multiplyBytes(dR, sA) + multiplyBytes(sR, BYTE_MAX_VALUE - dA);
    int resultG = multiplyBytes(dG, sA) + multiplyBytes(sG, BYTE_MAX_VALUE - dA);
    int resultB = multiplyBytes(dB, sA) + multiplyBytes(sB, BYTE_MAX_VALUE - dA);
    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}
static GPixel kXor(const GPixel source, const GPixel dest) {
    //!< [Sa + Da - 2 * Sa * Da, Sc * (1 - Da) + Dc * (1 - Sa)]
    int sA = GPixel_GetA(source);
    int sR = GPixel_GetR(source);
    int sG = GPixel_GetG(source);
    int sB = GPixel_GetB(source);
    
    int dA = GPixel_GetA(dest);
    int dR = GPixel_GetR(dest);
    int dG = GPixel_GetG(dest);
    int dB = GPixel_GetB(dest);

    int resultA = sA + dA - multiplyBytes(2 * sA, dA);
    int resultR = multiplyBytes(sR, BYTE_MAX_VALUE - dA) + multiplyBytes(dR, BYTE_MAX_VALUE - sA);
    int resultG = multiplyBytes(sG, BYTE_MAX_VALUE - dA) + multiplyBytes(dG, BYTE_MAX_VALUE - sA);
    int resultB = multiplyBytes(sB, BYTE_MAX_VALUE - dA) + multiplyBytes(dB, BYTE_MAX_VALUE - sA);
    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}
typedef GPixel(*ModeProc) (GPixel, GPixel);

class MyCanvas : public GCanvas {
public:
    MyCanvas(const GBitmap& device) : fDevice(device) {}
    
    /**
     *  Fill the entire canvas with the specified color, using the specified blendmode.
     */
    void drawPaint(const GPaint& paint) override {
        // convert paint input's color (GColor) to GPixel.
        const GColor& color = paint.getColor();
        GPixel srcP = convertColorToPixel(color);
        
        // get paint mode.
        const ModeProc procArray[12] = {kClear, kSrc, kDst, kSrcOver, kDstOver, kSrcIn,
                                    kDstIn, kSrcOut, kDstOut, kSrcATop, kDstATop, kXor};
        ModeProc proc = procArray[static_cast<int>(paint.getBlendMode())];

        // loop through the whole bitmap canvas, fill the color pixel by pixel.
        for (int y = 0; y < (this->fDevice).height(); ++y) {
            for (int x = 0; x < (this->fDevice).width(); ++x){
                GPixel* addr = (this->fDevice).getAddr(x, y);
                GPixel destP = *addr;
                *addr = proc(srcP, destP);
            }
        }
    }
    
    /**
     *  Fill the rectangle with the color, using the specified blendmode.
     *
     *  The affected pixels are those whose centers are "contained" inside the rectangle:
     *      e.g. contained == center > min_edge && center <= max_edge
     */
    void drawRect(const GRect& rect, const GPaint& paint) override {
        
        // clip: impose geometric restrictions
        GRect safeRect = GRect::MakeWH((this->fDevice).width(), (this->fDevice).height());
        if (!safeRect.intersect(rect)) {
            return;
        }
        
        // round the edges of the non-empty safeRect
        GIRect roundedSafeRect = safeRect.round();
        
        // get the source GColor from input paint and convert to Pixel
        const GColor& color = paint.getColor();
        GPixel srcP = convertColorToPixel(color);

        // get paint mode.
        const ModeProc procArray[12] = {kClear, kSrc, kDst, kSrcOver, kDstOver, kSrcIn,
                                    kDstIn, kSrcOut, kDstOut, kSrcATop, kDstATop, kXor};
        ModeProc proc = procArray[static_cast<int>(paint.getBlendMode())];

        for (int y = roundedSafeRect.top(); y < roundedSafeRect.bottom(); y++) {
            for (int x = roundedSafeRect.left(); x < roundedSafeRect.right(); x++) {
                GPixel* destAddr = (this->fDevice).getAddr(x, y);
                GPixel destP = *destAddr;
                GPixel resultP = proc(srcP, destP);
                *destAddr = resultP;
            }
        } 
    }

    /**
     *  Fill the convex polygon with the color, following the same "containment" rule as
     *  rectangles.
    */
    void drawConvexPolygon(const GPoint points[], int count, const GPaint& paint) override {
        // clip the points and make edges.
        if (count < 3) {
            return; // a polygon has to have at least 3 edges.
        }
        std::vector<Edge> edges;
        float bottom = (float)((this->fDevice).height());
        float right = (float)((this->fDevice).width());
        for (int i = 0; i < count - 1; ++i) {
            clipPolygonPoints(
                points[i], points[i+1], 
                0.0f, /* top */
                bottom,
                0.0f, /* left */
                right,
                edges);
        }
        clipPolygonPoints(points[count - 1], points[0], 0.0f, bottom, 0.0f, right, edges); // close the polygon.

        // terminate if no valid edges.
        if (edges.size() == 0) {
            return;
        }

        // sort by Ymin - Xmin - slope_min.
        std::sort(edges.begin(), edges.end());

        // get the source GColor from input paint and convert to Pixel
        const GColor& color = paint.getColor();
        GPixel srcP = convertColorToPixel(color);
        // get paint mode.
        const ModeProc procArray[12] = {kClear, kSrc, kDst, kSrcOver, kDstOver, kSrcIn,
                                    kDstIn, kSrcOut, kDstOut, kSrcATop, kDstATop, kXor};
        ModeProc proc = procArray[static_cast<int>(paint.getBlendMode())];

        // walk & blit.
        Edge leftE = edges[0];
        Edge rightE = edges[1];
        for (int y = edges[0].top_y, i = 1; y < edges.back().bottom_y; ++y) {
            if (y >= leftE.bottom_y) {
                leftE = edges[++i];
            }
            if (y >= rightE.bottom_y) {
                rightE = edges[++i];
            }
            for (int x = GRoundToInt(leftE.curr_x); x < GRoundToInt(rightE.curr_x); ++x) {
                GPixel* destAddr = (this->fDevice).getAddr(x, y);
                GPixel destP = *destAddr;
                GPixel resultP = proc(srcP, destP);
                *destAddr = resultP;
            }
            leftE.curr_x = leftE.slope + leftE.curr_x;
            rightE.curr_x = rightE.slope + rightE.curr_x;
        }
    }
    
private:
    const GBitmap fDevice;
};

std::unique_ptr<GCanvas> GCreateCanvas(const GBitmap& device) {
    if (!device.pixels()) {
        return nullptr;
    }
    return std::unique_ptr<GCanvas>(new MyCanvas(device));
}

