//
// template-rt.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

#define MAX_DISTANCE 10000.0f
#define MIN_DIST 1.0f
#define MIN_REFLECT_DIST 0.001f
#define MAX_REFLECTION_DEPTH 3


struct Ray
{
    vec4 origin;
    vec4 dir;
    bool hasReflected;
};

struct Sphere
{
    string  name;
    vec4    position;
    vec3    scale;
    vec4    color;
    float   Ka,
    Kd,
    Ks,
    Kr,
    specularExponent;
    mat4    transform,
    inverseTransform;
    
    Sphere(string tname, vec4 tposition, vec3 tscale, vec4 tcolor, float tKa, float tKd, float tKs, float tKr, float tspe)
    {
        name = tname;
        position = tposition;
        scale = tscale;
        color = tcolor;
        Ka = tKa;
        Kd = tKd;
        Ks = tKs;
        Kr = tKr;
        specularExponent = tspe;
        transform = Translate(position) * Scale(scale);
        InvertMatrix(Scale(scale), inverseTransform);
    }
};

struct Light
{
    vec4 color;
    string name;
    vec4 position;
    
    Light(string _name, vec4 _position, vec4 _color)
    {
        name = _name;
        position = _position;
        color = _color;
    }
};

struct Intersection
{
    Ray ray;
    float distance;
    vec4 point;
    bool interiorPoint;
    Sphere *sphere;
    vec4 normal;
};

vector<vec4> g_colors;

int g_width;
int g_height;

float g_left;
float g_right;
float g_top;
float g_bottom;
float g_near;

std::vector<Sphere> g_spheres;
std::vector<Light> g_lights;

vec4 g_backgroundColor;
vec4 g_ambient;
char* g_outPutFileName;



// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}

void parseLine(const vector<string>& vs)
{
    // add parsing of NEAR, LEFT, RIGHT, BOTTOM, TOP, SPHERE, LIGHT, BACK, AMBIENT, OUTPUT.
    if(vs[0] == "NEAR")
        g_near = toFloat(vs[1]);
    else if(vs[0] == "LEFT")
        g_left = toFloat(vs[1]);
    else if(vs[0] == "RIGHT")
        g_right = toFloat(vs[1]);
    else if(vs[0] == "BOTTOM")
        g_bottom = toFloat(vs[1]);
    else if(vs[0] == "TOP")
        g_top = toFloat(vs[1]);
    else if (vs[0] == "RES")
    {
        g_width = (int)toFloat(vs[1]);
        g_height = (int)toFloat(vs[2]);
        g_colors.resize(g_width * g_height);
    }
    else if(vs[0] == "SPHERE")
    {
        if(g_spheres.size() < 5)
            g_spheres.push_back(Sphere(vs[1], // name
                                       toVec4(vs[2], vs[3], vs[4]), // position
                                       vec3(toFloat(vs[5]), toFloat(vs[6]), toFloat(vs[7])), // scale
                                       toVec4(vs[8], vs[9], vs[10]), // color
                                       toFloat(vs[11]), // Ka
                                       toFloat(vs[12]), // Kd
                                       toFloat(vs[13]), // Ks
                                       toFloat(vs[14]), // Kr
                                       toFloat(vs[15]))); // specularExponent
    }
    else if(vs[0] == "LIGHT")
    {
        if(g_lights.size() < 5)
            g_lights.push_back(Light(vs[1], // name
                                     toVec4(vs[2], vs[3], vs[4]), // position
                                     toVec4(vs[5], vs[6], vs[7]))); // color
    }
    else if(vs[0] == "BACK")
        g_backgroundColor = toVec4(vs[1], vs[2], vs[3]);
    else if(vs[0] == "AMBIENT")
        g_ambient = toVec4(vs[1], vs[2], vs[3]);
    else if(vs[0] == "OUTPUT")
    {
        std::string str = vs[1];
        g_outPutFileName = new char[str.size() + 1];
        std::copy(str.begin(), str.end(), g_outPutFileName);
        g_outPutFileName[str.size()] = '\0';
    }
    else {}
}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
    /*
     cout << "g_width " << g_width << endl;
     cout << "g_height " << g_height << endl;
     cout << "g_left " << g_left << endl;
     cout << "g_right " << g_right << endl;
     cout << "g_top " << g_top << endl;
     cout << "g_bottom " << g_bottom << endl;
     cout << "g_near " << g_near << endl;
     */
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}

vec3 toVec3(vec4 in)
{
    return vec3( in[0], in[1], in[2] );
}


// -------------------------------------------------------------------
// Intersection routine

// TODO: add your ray-sphere intersection routine here.
Intersection getIntersection(Ray& ray)
{
    Intersection intersection;
    intersection.ray = ray;
    intersection.distance = -1;
    intersection.interiorPoint = false;
    
    for (Sphere &sphere : g_spheres) {
        vec4 S = sphere.inverseTransform * (sphere.position - ray.origin); // -(O - C)
        vec4 C = sphere.inverseTransform * ray.dir;
        
        // Quadratic equation: |c|^2t^2 + 2(S.tc) + |S|^2 - 1
        float a = dot(C, C);
        float b = dot(S, C);
        float c = dot(S, S) - 1;
        
        // Solve equation
        float solution;
        float discriminant = b * b - a * c; // Value under the root
        
        bool interiorPoint = false;
        
        if (discriminant < 0) {
            // No solutions: line does not intersect
            continue;
        } else if (discriminant == 0) {
            // Single solution: line intersects at one point
            solution = b / a;
        } else {
            // Two solutions: line intersects at two points
            float root = sqrtf(discriminant);
            float solution1 = (b - root) / a;
            float solution2 = (b + root) / a;
            
            // Use the smallest valid solution
            solution = fminf(solution1, solution2);
            if (solution <= MIN_REFLECT_DIST || ((!ray.hasReflected) && solution <= MIN_DIST)) {
                solution = fmaxf(solution1, solution2);
                interiorPoint = true;
            }
        }
        
        // Validate solution
        if (solution <= MIN_REFLECT_DIST || ((!ray.hasReflected) && solution <= MIN_DIST)) {
            continue;
        }
        
        if ((intersection.distance == -1 || solution < intersection.distance)) {
            intersection.distance = solution;
            intersection.sphere = &sphere;
            intersection.interiorPoint = interiorPoint;
        }
    }
    
    // Calculate intersection point and normal
    if (intersection.distance != -1) {
        ray.hasReflected = true;
        intersection.point = ray.origin + ray.dir * intersection.distance;
        
        vec4 normal = intersection.point - intersection.sphere->position;
        // Invert normal for interior points
        if (intersection.interiorPoint) {
            normal = -normal;
        }
        
        mat4 trans = transpose(intersection.sphere->inverseTransform);
        normal = trans * intersection.sphere->inverseTransform * normal;
        normal.w = 0;
        intersection.normal = normalize(normal);
    }
    
    return intersection;
}


// -------------------------------------------------------------------
// Ray tracing

vec4 checkColor(vec4 color)
{
    vec4 result = color;
    if (color.x < 0.0f) { result.x = 0.0f; }
    if (color.y < 0.0f) { result.y = 0.0f; }
    if (color.z < 0.0f) { result.z = 0.0f; }
    if (color.x > 1.0f) { result.x = 1.0f; }
    if (color.y > 1.0f) { result.y = 1.0f; }
    if (color.z > 1.0f) { result.z = 1.0f; }
    
    return result;
}



vec4 trace(Ray& ray, int reflectionDepth)
{
    // Limit reflection level
    if (reflectionDepth >= MAX_REFLECTION_DEPTH) {
        return vec4();
    }
    
    // Get the nearest intersecting sphere
    Intersection intersection = getIntersection(ray);
    
    if (intersection.distance == -1 && reflectionDepth == 0) {
        // Return background color if no intersection and is an initial ray
        return g_backgroundColor;
    } else if (intersection.distance == -1) {
        // Return no color if no intersection and not an initial ray
        return vec4();
    }
    
    // Calculate initial intersection color with ambient intensity
    vec4 color = intersection.sphere->color * intersection.sphere->Ka * g_ambient;
    
    // Calculate Blinn-Phong shading from light sources
    vec4 diffusion = vec4(0, 0, 0, 0);
    vec4 specular = vec4(0, 0, 0, 0);
    for (Light light : g_lights) {
        // Generate ray from intersection point to light
        Ray lightRay;
        lightRay.origin = intersection.point;
        lightRay.dir = normalize(light.position - intersection.point);
        
        // Determine if the light source is not obstructed
        Intersection lightIntersection = getIntersection(lightRay);
        if (lightIntersection.distance == -1) {
            // Calculate the intensity of diffuse light
            float diffusionIntensity = dot(intersection.normal, lightRay.dir);
            if (diffusionIntensity > 0) {
                diffusion += diffusionIntensity * light.color * intersection.sphere->color;
                
                // Calculate the half vector between light vector and the view vector
                vec4 H = normalize(lightRay.dir - ray.dir);
                
                // Calculate the intensity of specular light
                float specularIntensity = dot(intersection.normal, H);
                specular += powf(powf(specularIntensity, intersection.sphere->specularExponent), 3) * light.color;
            }
        }
    }
    
    // Apply diffusion and specular values
    color += diffusion * intersection.sphere->Kd + specular * intersection.sphere->Ks;
    
    // Calculate reflections
    Ray reflectRay;
    reflectRay.origin = intersection.point;
    reflectRay.dir = normalize(ray.dir - 2.0f * intersection.normal * dot(intersection.normal, ray.dir));
    reflectionDepth += 1;

    color += trace(reflectRay, reflectionDepth) * intersection.sphere->Kr;
    
    return color;
}


vec4 getDir(int ix, int iy)
{
    // This should return the direction from the origin
    // to pixel (ix, iy), normalized.
    
    float alpha = (float) ix / g_width;
    float beta = (float) iy / g_height;
    
    float x = (1.0f - alpha) * g_left + alpha * g_right;
    float y = (1.0f - beta) * g_bottom + beta * g_top;
    
    //cout << "x " << x << "   y " << y << endl;
    
    vec4 dir;
    dir = vec4(x, y, -g_near, 0.0f);
    return dir;
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    ray.hasReflected = false;
    int reflectionDepth = 0;
    vec4 color = trace(ray, reflectionDepth);
    setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, char* fname, unsigned char* pixels)
{
    FILE *fp;
    const int maxVal=255;
    
    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);
    
    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }
    
    fclose(fp);
}

float clamp(float val) {
    if(val > 1.0f)
        return 1.0f;
    else if( val < 0.0f)
        return 0.0f;
    else
        return val;
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    // Clamp values if out of range.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++)
            {
                float clampedValue = clamp(((float*)g_colors[y*g_width+x])[i]);
                buf[y*g_width*3+x*3+i] = (unsigned char)(clampedValue * 255.9f);
            }
    
    // Change file name based on input file name.
    savePPM(g_width, g_height, g_outPutFileName, buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
    render();
    saveFile();
    return 0;
}
