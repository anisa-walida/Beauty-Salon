#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#define STB_IMAGE_IMPLEMENTATION
#pragma warning(disable:4996)
#include <unordered_map>
#include "shader.h"
#include "camera.h"
#include "basic_camera.h"
#include "pointLight.h"
#include "cube.h"
#include "Cone.h"
#include "sphere.h"
#include "cylinder.h"
#include "stb_image.h"
#include <iostream>
#include <ctime>


using namespace std;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);
void drawCube(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 model, float r, float g, float b);
void getCurrentTime(int& hours, int& minutes, int& seconds);
glm::mat4 transform(float tx, float ty, float tz, float sx, float sy, float sz);
// ************************Waiting ROOM****************************************
void f_waiting_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube);
void f_waiting_sofa(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& sofa_top, Cube& sofa_foam, Cube& sofa_pillow);
void f_waiting_sofa1(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& sofa_top, Cube& sofa_foam, Cube& sofa_pillow, float rotateY, float tx, float ty, float tz);
void f_waiting_table(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_table);
void f_waiting_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& wall, Cube& door,Cube& sunroof, Cube& glass, Cube& green, Shader& ourShader);
void f_waiting_tv(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_tv, Cube& drawing_sound_box, Cube& drawing_cupboard);
void f_waiting_mat(unsigned int& cylinderVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cylinder& drawing_mat);
void fan(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& b, Cube& c, float x = 0.0f, float y = 0.0f, float z = 0.0f);
void f_waiting_window(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_window);


//  *************************************ROOM 1******************************************
void f_room1_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube, Shader& ourShader, Sphere& clock_bell, Cube& clock, Cube& x);
void f_room1_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube);
void f_room1_bed(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& bed_sheet, Cube& bed_pillow, Cube& bed_texture, Cube& blanket_texture);
void f_room1_almari(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& draw_almari,Cube& door_right,Cube& door_left);
void f_room1_window(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_bed_room1_window);
void f_room1_bedtable(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_bed_room1_bedtable, Shader& ourShader);
//  *************************************ROOM2******************************************
void f_room2_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube);
void f_room2_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube);
void f_room2_table(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& dining_table, Cube& mirror);
void f_room2_table2(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& dining_table, Cube& mirror);
void f_room2_frize(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_dining_frize);
void f_room2_window(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& kitchen_window);
void f_room2_light(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Shader& ourShader);
void chair(unsigned int& cubeVAO, Shader& lightingShaderWithTexture, glm::mat4 alTogether, Cube& chair_texture,Cube& dining_frize,float tx, float ty, float tz, float rotateY);
void chair2(unsigned int& cubeVAO, Shader& lightingShaderWithTexture, glm::mat4 alTogether, Cube& chair_texture, Cube& dining_frize, float tx, float ty, float tz, float rotateY);


//  *************************************ROOM 3********************************************
void f_room3_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube);
void f_room3_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube);
void f_room3_book(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_dining_frize);
void f_room3_light(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Shader& ourShader);
void f_room3_dressing_table(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& dressing_right, Cube& dressing_bottom, Cube& dressing_mirror);
void f_room3_bed(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& bed_sheet, Cube& bed_pillow, Cube& bed_texture, Cube& blanket_texture);

// **************************************ROOF***************************************************
void f_roof(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_roof);
unsigned int loadTexture(char const* path, GLenum textureWrappingModeS, GLenum textureWrappingModeT, GLenum textureFilteringModeMin, GLenum textureFilteringModeMax);
void shaderActivate(Shader& shader);
void f_bathroom_bucket(unsigned int& cubeVAO,Shader& lightingShader,glm::mat4 alTogether,Shader& lightingShaderWithTexture,Cube& water);



//////////Plant
void drawFractalCactus(unsigned int& cubeVAO, Shader& lightingShaderWithTexture, glm::mat4 alTogether,
    Cube& base_texture, Cube& cactus_texture,
    float tx, float ty, float tz, int depth, float scaleFactor, float branchAngle);


//////////Pedestal
void drawPedestal(Shader& lightingShaderWithTexture
    ,unsigned int diffuseMap, unsigned int specularMap);

/////////Chandelier
void chandelier(
    unsigned int& cubeVAO,Shader& lightingShader,glm::mat4 alTogether,Shader& lightingShaderWithTexture,Cube& cube,Sphere& sphere,unsigned int diffMap8,unsigned int specMap8,float x,float y,float z);


void drawMirrorWithStand(Shader& lightingShaderWithTexture,glm::mat4 alTogether,Cube& frame,Cube& mirrorSurface, Cube& stand,float tx,float ty,float tz,float rotateY,float scaleFactor);

void duvan(Cube& cube2, Cube& cube3, Shader& lightingShaderWithTexture, glm::mat4 als, unsigned int diffuseMap, unsigned int specularMap);
void drawLipstick(Shader& shader, glm::mat4 baseTransform, unsigned int diffuseMap, unsigned int specularMap, unsigned int diffuseMapcone, unsigned int specularMapcone );

void drawShowpiece(Shader& shader, glm::mat4 baseTransform, unsigned int diffuseMap, unsigned int specularMap, unsigned int diffuseMapcone, unsigned int specularMapcone, float rotationAngle) {
    glm::mat4 model;
    glm::mat4 translateMatrix, scaleMatrix, rotateMatrix;

    // Apply uniform scaling and translation to the entire showpiece
    float uniformScale = 0.5f; // Scale the showpiece to half its size
    glm::vec3 uniformTranslation = glm::vec3(1.5f, 0.2f, -1.0f); // Translate the entire showpiece to a new position
    baseTransform = glm::translate(baseTransform, uniformTranslation);
    baseTransform = glm::scale(baseTransform, glm::vec3(uniformScale, uniformScale, uniformScale));

    // Apply rotation to the entire showpiece
    // This makes the entire showpiece rotate around the center of the base cylinder
    glm::vec3 rotationCenter = glm::vec3(3.0f, 0.7f, -2.0f + 0.2f); // Center of rotation (center of base cylinder)
    baseTransform = glm::translate(baseTransform, rotationCenter); // Move to center
    baseTransform = glm::rotate(baseTransform, glm::radians(rotationAngle), glm::vec3(0.0f, 1.0f, 0.0f)); // Rotate around Y-axis
    baseTransform = glm::translate(baseTransform, -rotationCenter); // Move back to original position

    // Create the base cylinder
    Cylinder baseCylinder(0.5f, 0.2f, 36);  // Base cylinder

    // Separate transformation for the base cylinder
    glm::mat4 baseCylinderTransform = glm::translate(glm::mat4(1.0f), glm::vec3(3.0f, 0.7f, -2.0f + 0.2f)); // Base is slightly below lipsticks
    baseCylinderTransform = glm::scale(baseCylinderTransform, glm::vec3(0.7f, 1.0f, 0.7f)); // Scale base cylinder slightly larger
    model = baseTransform * baseCylinderTransform; // Combine baseTransform with the cylinder's specific transform
    shader.setMat4("model", model);
    baseCylinder.drawCylinderTexture(shader, model, diffuseMap, specularMap);

    // Predefined positions for the lipsticks around the cylinder
    std::vector<glm::vec3> lipstickPositions = {
        glm::vec3(0.3f, 0.5f, 0.0f),   // Lipstick 1
        glm::vec3(0.1f, 0.5f, 0.25f), // Lipstick 2
        glm::vec3(-0.2f, 0.5f, 0.15f), // Lipstick 3
        glm::vec3(-0.2f, 0.5f, -0.15f), // Lipstick 4
        glm::vec3(0.1f, 0.5f, -0.25f)  // Lipstick 5
    };

    // Draw each lipstick at its predefined position
    for (int i = 0; i < lipstickPositions.size(); ++i) {
        // Translate lipstick to its position
        translateMatrix = glm::translate(glm::mat4(1.0f), lipstickPositions[i]);

        // Combine transformations
        glm::mat4 lipstickTransform = baseTransform * translateMatrix;

        // Draw lipstick
        drawLipstick(shader, lipstickTransform, diffuseMap, specularMap, diffuseMapcone, specularMapcone);

        // Print positions for debugging
        glm::vec4 position = baseTransform * translateMatrix * glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
        // std::cout << "Lipstick position " << (i + 1) << ": " << position.x << ", " << position.y << ", " << position.z << std::endl;
    }
}



// settings
const unsigned int SCR_WIDTH = 1000;
const unsigned int SCR_HEIGHT = 800;

// modelling transform
float rotateAngle_X = 0.0;
float rotateAngle_Y = 0.0;
float rotateAngle_Z = 0.0;
float rotateAxis_X = 0.0;
float rotateAxis_Y = 0.0;
float rotateAxis_Z = 1.0;
float translate_X = 0.0;
float translate_Y = 0.0;
float translate_Z = 0.0;
float scale_X = 1.0;
float scale_Y = 1.0;
float scale_Z = 1.0;

// front door
float f_door = 0.0f;
float b1_door = 0.0f;
// bathroom door
float max_bathroom_door_translate = 0.6;
float min_bathroom_door_translate = -0.6;
float bathroom_door_translate = 0.0f;

// fan
bool isFanOn = false;
float rotateFan = 0;
float rotateClock = 0.0f;
bool sign = 1;

//chair rotation
float chairRotation = 0.0f;
float backrestTiltAngle = -10.0f;
//sponge
float spongeup = 0.0f;
//showpiece
float rotationAngle = 0.0f; // Initialize rotation angle
float rotationSpeed = 50.0f; // Rotation speed (degrees per second)


// camera
Camera camera(glm::vec3(0.0f, 1.1f, 5.2f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

float eyeX = 0.0, eyeY = 1.0, eyeZ = 3.0;
float lookAtX = 0.0, lookAtY = 0.0, lookAtZ = 0.0;
glm::vec3 V = glm::vec3(0.0f, 1.0f, 0.0f);
BasicCamera basic_camera(eyeX, eyeY, eyeZ, lookAtX, lookAtY, lookAtZ, V);



// positions of the point lights
glm::vec3 pointLightPositions[] = {
    glm::vec3(0.17f, 0.4f, -1.75f),
    glm::vec3(0.0f,  1.5f,  0.0f),
    glm::vec3(0.0f,  1000.0f,  0.0f),
    glm::vec3(0.0f,  3.0f,  0.0f)
};

glm::vec3 point_light_positions[] = {
    glm::vec3(1.45f, 1.3f, 0.1f),
    glm::vec3(1.45f, 1.3f, -3.1f),
    glm::vec3(1.6f, 1.3f, -3.1f),
    glm::vec3(1.6f, 1.3f, 0.5f),
    glm::vec3(2.5f + 1.9f, 0.8f, -0.9f)
};

PointLight pointlight1(

    pointLightPositions[0].x, pointLightPositions[0].y, pointLightPositions[0].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    0.7f, 0.7f, 0.7f,        // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    1       // light number
);
PointLight pointlight2(

    pointLightPositions[1].x, pointLightPositions[1].y, pointLightPositions[1].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    0.7f, 0.7f, 0.7f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    2       // light number
);

PointLight pointlight3(

    pointLightPositions[2].x, pointLightPositions[2].y, pointLightPositions[2].z,  // position
    0.1f, 0.1f, 0.1f,     // ambient
    0.1f, 0.1f, 0.1f,      // diffuse
    0.1f, 0.1f, 0.1f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    3       // light number
);
PointLight pointlight4(

    pointLightPositions[3].x, pointLightPositions[3].y, pointLightPositions[3].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    0.7f, 0.7f, 0.7f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    4       // light number
);
// ******************************DRAWING_ROOM_LIGHT***********************************
PointLight drawing_light(
    point_light_positions[0].x, point_light_positions[0].y, point_light_positions[0].z,  // position
    0.3f, 0.3f, 0.3f,     // ambient
    0.3f, 0.3f, 0.3f,      // diffuse
    0.3f, 0.3f, 0.3f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    5       // light number
);
PointLight bed_room1_light(
    point_light_positions[1].x, point_light_positions[1].y, point_light_positions[1].z,  // position
    0.3f, 0.3f, 0.3f,     // ambient
    0.3f, 0.3f, 0.3f,      // diffuse
    0.3f, 0.3f, 0.3f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    6       // light number
);
PointLight dining_light(
    point_light_positions[2].x, point_light_positions[2].y, point_light_positions[2].z,  // position
    0.1f, 0.1f, 0.1f,     // ambient
    0.2f, 0.2f, 0.2f,      // diffuse
    0.2f, 0.2f, 0.2f,         // specular
    1.0f,   //k_c
    0.14f,  //k_l
    0.07f,  //k_q
    7       // light number
);
PointLight bed_room2_light(
    point_light_positions[3].x, point_light_positions[3].y, point_light_positions[3].z,  // position
    0.3f, 0.3f, 0.3f,     // ambient
    0.3f, 0.3f, 0.3f,      // diffuse
    0.3f, 0.3f, 0.3f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    8       // light number
);
PointLight bathroom_light(
    point_light_positions[4].x, point_light_positions[4].y, point_light_positions[4].z,  // position
    0.2f, 0.2f, 0.2f,     // ambient
    0.2f, 0.2f, 0.2f,      // diffuse
    0.2f, 0.2f, 0.2f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    9       // light number
);



// light settings
bool onOffPointToggle = true;
bool onOffSpotToggle = true;
bool onOffDirectToggle = true;
bool ambientToggle = true;
bool diffuseToggle = true;
bool specularToggle = true;

//glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
//glm::mat4 view = camera.GetViewMatrix();
glm::mat4 projection;
glm::mat4 view;

string diffuseMapPath;
string specularMapPath;

class Curve
{
public:
    vector<float> cntrlPoints;
    vector <float> coordinates;
    vector <float> normals;
    vector <int> indices;
    vector <float> vertices;
    const double pi = 3.14159265389;
    const int nt = 40;
    const int ntheta = 20;
    Curve(vector<float>& tmp)
    {
        this->cntrlPoints = tmp;
        this->fishVAO = hollowBezier(cntrlPoints.data(), ((unsigned int)cntrlPoints.size() / 3) - 1);
        cout << cntrlPoints.size() << endl;
        cout << coordinates.size() << endl;
        cout << normals.size() << endl;
        cout << indices.size() << endl;
        cout << vertices.size() << endl;
    }
    ~Curve()
    {
        glDeleteVertexArrays(1, &fishVAO);
        glDeleteVertexArrays(1, &bezierVAO);
        glDeleteBuffers(1, &bezierVBO);
        glDeleteBuffers(1, &bezierEBO);
    }
    void draw(Shader& lightingShader, glm::mat4 model)
    {
        /// Fish
        lightingShader.use();
        lightingShader.setMat4("model", model);
        lightingShader.setVec3("material.ambient", glm::vec3(1.0f, 0.6f, 0.0f));
        lightingShader.setVec3("material.diffuse", glm::vec3(1.0f, 0.6f, 0.0f));
        lightingShader.setVec3("material.specular", glm::vec3(1.0f, 1.0f, 1.0f));
        lightingShader.setFloat("material.shininess", 32.0f);

        glBindVertexArray(fishVAO);
        glDrawElements(GL_TRIANGLES,                    // primitive type
            (unsigned int)indices.size(),          // # of indices
            GL_UNSIGNED_INT,                 // data type
            (void*)0);                       // offset to indices

        // unbind VAO
        glBindVertexArray(0);
        /// End Fish
    }
private:
    unsigned int fishVAO;
    unsigned int bezierVAO;
    unsigned int bezierVBO;
    unsigned int bezierEBO;


    unsigned int drawControlPoints()
    {
        unsigned int controlPointVAO;
        unsigned int controlPointVBO;

        glGenVertexArrays(1, &controlPointVAO);
        glGenBuffers(1, &controlPointVBO);

        glBindVertexArray(controlPointVAO);

        glBindBuffer(GL_ARRAY_BUFFER, controlPointVBO);
        glBufferData(GL_ARRAY_BUFFER, (unsigned int)cntrlPoints.size() * sizeof(float), cntrlPoints.data(), GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        return controlPointVAO;
    }

    long long nCr(int n, int r)
    {
        if (r > n / 2)
            r = n - r; // because C(n, r) == C(n, n - r)
        long long ans = 1;
        int i;

        for (i = 1; i <= r; i++)
        {
            ans *= n - r + i;
            ans /= i;
        }

        return ans;
    }
    void BezierCurve(double t, float xy[2], GLfloat ctrlpoints[], int L)
    {
        double y = 0;
        double x = 0;
        t = t > 1.0 ? 1.0 : t;
        for (int i = 0; i < L + 1; i++)
        {
            long long ncr = nCr(L, i);
            double oneMinusTpow = pow(1 - t, double(L - i));
            double tPow = pow(t, double(i));
            double coef = oneMinusTpow * tPow * ncr;
            x += coef * ctrlpoints[i * 3];
            y += coef * ctrlpoints[(i * 3) + 1];

        }
        xy[0] = float(x);
        xy[1] = float(y);
    }
    unsigned int hollowBezier(GLfloat ctrlpoints[], int L)
    {
        int i, j;
        float x, y, z, r;                //current coordinates
        float theta;
        float nx, ny, nz, lengthInv;    // vertex normal


        const float dtheta = 2 * pi / ntheta;        //angular step size

        float t = 0;
        float dt = 1.0 / nt;
        float xy[2];

        for (i = 0; i <= nt; ++i)              //step through y
        {
            BezierCurve(t, xy, ctrlpoints, L);
            r = xy[0];
            y = xy[1];
            theta = 0;
            t += dt;
            lengthInv = 1.0 / r;

            for (j = 0; j <= ntheta; ++j)
            {
                double cosa = cos(theta);
                double sina = sin(theta);
                z = r * cosa;
                x = r * sina;

                coordinates.push_back(x);
                coordinates.push_back(y);
                coordinates.push_back(z);

                // normalized vertex normal (nx, ny, nz)
                // center point of the circle (0,y,0)
                nx = (x - 0) * lengthInv;
                ny = (y - y) * lengthInv;
                nz = (z - 0) * lengthInv;

                normals.push_back(nx);
                normals.push_back(ny);
                normals.push_back(nz);

                theta += dtheta;
            }
        }
        // generate index list of triangles
        // k1--k1+1
        // |  / |
        // | /  |
        // k2--k2+1

        int k1, k2;
        for (int i = 0; i < nt; ++i)
        {
            k1 = i * (ntheta + 1);     // beginning of current stack
            k2 = k1 + ntheta + 1;      // beginning of next stack

            for (int j = 0; j < ntheta; ++j, ++k1, ++k2)
            {
                // k1 => k2 => k1+1
                indices.push_back(k1);
                indices.push_back(k2);
                indices.push_back(k1 + 1);

                // k1+1 => k2 => k2+1
                indices.push_back(k1 + 1);
                indices.push_back(k2);
                indices.push_back(k2 + 1);
            }
        }

        size_t count = coordinates.size();
        for (int i = 0; i < count; i += 3)
        {
            //cout << count << ' ' << i + 2 << endl;
            vertices.push_back(coordinates[i]);
            vertices.push_back(coordinates[i + 1]);
            vertices.push_back(coordinates[i + 2]);

            vertices.push_back(normals[i]);
            vertices.push_back(normals[i + 1]);
            vertices.push_back(normals[i + 2]);
        }

        glGenVertexArrays(1, &bezierVAO);
        glBindVertexArray(bezierVAO);

        // create VBO to copy vertex data to VBO
        glGenBuffers(1, &bezierVBO);
        glBindBuffer(GL_ARRAY_BUFFER, bezierVBO);           // for vertex data
        glBufferData(GL_ARRAY_BUFFER,                   // target
            (unsigned int)vertices.size() * sizeof(float), // data size, # of bytes
            vertices.data(),   // ptr to vertex data
            GL_STATIC_DRAW);                   // usage

        // create EBO to copy index data
        glGenBuffers(1, &bezierEBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bezierEBO);   // for index data
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,           // target
            (unsigned int)indices.size() * sizeof(unsigned int),             // data size, # of bytes
            indices.data(),               // ptr to index data
            GL_STATIC_DRAW);                   // usage

        // activate attrib arrays
        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);

        // set attrib arrays with stride and offset
        int stride = 24;     // should be 24 bytes
        glVertexAttribPointer(0, 3, GL_FLOAT, false, stride, (void*)0);
        glVertexAttribPointer(1, 3, GL_FLOAT, false, stride, (void*)(sizeof(float) * 3));

        // unbind VAO, VBO and EBO
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

        return bezierVAO;
    }

};

Curve* bucket;
Curve* backrest;
Curve* seat;
Curve* backrest2;

vector<float> Bucket = {
  

    // Front top edge (backrest curve)
 -1.0, 1.0, 0.2,  // Left
 -0.5, 1.5, 0.2,  // Curve point
  0.5, 1.5, 0.2,  // Curve point
  1.0, 1.0, 0.2,  // Right
  // Front top edge (backrest curve)
  1.0, 1.0, 0.2,  // Left
  0.5, 1.5, 0.2,  // Curve point
  -0.5, 1.5, 0.2,  // Curve point
  -1.0, 1.0, 0.2,  // Right
     
};
vector<float> Backrest = {


    // Front top edge (backrest curve)
  -1.0, 1.0, 0.2,  // Left
  -0.5, 1.5, 0.2,  // Curve point
   0.5, 1.5, 0.2,  // Curve point
   1.0, 1.0, 0.2,  // Right
   // Front top edge (backrest curve)
   1.0, 1.0, 0.2,  // Left
   0.5, 1.5, 0.2,  // Curve point
   -0.5, 1.5, 0.2,  // Curve point
   -1.0, 1.0, 0.2,  // Right

};
vector<float> Seat = {
    // Front Face (Curved)
    -0.65, 0.0,  0.2,  // Left-bottom
    -0.5, 1.5,  0.4,  // Left-curve
     0.5, 1.5,  0.4,  // Right-curve
     0.65, 0.0,  0.2,  // Right-bottom

     // Back Face (Curved)
     -0.65, 0.0, -0.2,  // Left-bottom
     -0.5, -1.5, -0.2,  // Left-curve
      0.5, -1.5, -0.2,  // Right-curve
      0.65, 0.0, -0.2,  // Right-bottom

};
vector<float> Backrest2 = {
    // Front Face (Curved)
    -0.65, 0.0,  0.2,  // Left-bottom
    -0.5, 1.5,  0.4,  // Left-curve
     0.5, 1.5,  0.4,  // Right-curve
     0.65, 0.0,  0.2,  // Right-bottom

     // Back Face (Curved)
     -0.65, 0.0, -0.2,  // Left-bottom
     -0.5, -1.5, -0.2,  // Left-curve
      0.5, -1.5, -0.2,  // Right-curve
      0.65, 0.0, -0.2,  // Right-bottom

};


// timing
float deltaTime = 0.0f;    // time between current frame and last frame
float lastFrame = 0.0f;

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Beauty Salon", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // tell GLFW to capture our mouse
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);

    // build and compile our shader zprogram
    // ------------------------------------
    Shader lightingShader("vertexShaderForPhongShading.vs", "fragmentShaderForPhongShading.fs");
    Shader lightingShaderWithTexture("vertexShaderForPhongShadingWithTexture.vs", "fragmentShaderForPhongShadingWithTexture.fs");
    Shader ourShader("vertexShader.vs", "fragmentShader.fs");

    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------

    float cube_vertices[] = {
        // positions      // normals
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        1.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        1.0f, 1.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        0.0f, 1.0f, 0.0f, 0.0f, 0.0f, -1.0f,

        1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,

        0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
        1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
        1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f,
        0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f,

        0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 1.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f,

        1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f,
        1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f,

        0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
        1.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f
    };
    unsigned int cube_indices[] = {
        0, 3, 2,
        2, 1, 0,

        4, 5, 7,
        7, 6, 4,

        8, 9, 10,
        10, 11, 8,

        12, 13, 14,
        14, 15, 12,

        16, 17, 18,
        18, 19, 16,

        20, 21, 22,
        22, 23, 20
    };

    unsigned int cubeVAO, cubeVBO, cubeEBO,cylinderVAO,sphereVAO;
    glGenVertexArrays(1, &cubeVAO);
    glGenBuffers(1, &cubeVBO);
    glGenBuffers(1, &cubeEBO);

    glBindVertexArray(cubeVAO);

    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cube_vertices), cube_vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cube_indices), cube_indices, GL_STATIC_DRAW);


    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // vertex normal attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)12);
    glEnableVertexAttribArray(1);

    // second, configure the light's VAO (VBO stays the same; the vertices are the same for the light object which is also a 3D cube)
    unsigned int lightCubeVAO;
    glGenVertexArrays(1, &lightCubeVAO);
    glBindVertexArray(lightCubeVAO);

    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeEBO);
    // note that we update the lamp's position attribute's stride to reflect the updated buffer data
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);


    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // drawing_floor--------------------------------------------------------

    // Declare a Sphere instance globally or within the main function
    // Define texture paths
    

    // Load textures
    unsigned int diffuseMap = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specularMap = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);

    

    diffuseMapPath = "images/floor3.jpg";
    specularMapPath = "images/floor3.jpg";
    Cube drawing_floor = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 6.0f, 4.0f);

    // fan_center-----------------------------------------------------------
    diffuseMapPath = "images/fan_center.PNG";
    specularMapPath = "images/fan_center.PNG";
    Cube c = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // fan_blade--------------------------------------------------------------------
    diffuseMapPath = "images/f_b.png";
    specularMapPath = "images/f_b.png";
    Cube b = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // drawing_wall-----------------------------------------------------------
    diffuseMapPath = "images/wall_color.jpg";
    specularMapPath = "images/wall_color.jpg";
    Cube drawing_wall = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 6.0f, 4.0f);

    // drawing_window-----------------------------------------------------------
    diffuseMapPath = "images/drawing_window.jpg";
    specularMapPath = "images/drawing_window.jpg";
    Cube drawing_window = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    //   sofa_top ------------------------------------------------------------
    diffuseMapPath = "images/sofa_top.jpg";
    specularMapPath = "images/sofa_top.jpg";
    Cube sofa_top = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 6.0f, 4.0f);


    // sofa_foam --------------------------------------------------------------
    diffuseMapPath = "images/sofa_foam.jpg";
    specularMapPath = "images/sofa_foam.jpg";
    Cube sofa_foam = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 6.0f, 4.0f);

    // sofa_pillow-----------------------------------------------------------------
    diffuseMapPath = "images/sofa_pillow.jpg";
    specularMapPath = "images/sofa_pillow.jpg";
    Cube sofa_pillow = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 6.0f, 4.0f);

    //  drawing_table---------------------------------------------------------------
    diffuseMapPath = "images/table_top.jpg";
    specularMapPath = "images/table_top.jpg";
    Cube drawing_table = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 6.0f, 4.0f);

    //  drawing_tv-------------------------------------------------------------------
    diffuseMapPath = "images/tv.jpg";
    specularMapPath = "images/tv.jpg";
    Cube drawing_tv = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    //  drawing_sound_box-------------------------------------------------------------
    diffuseMapPath = "images/soundbox.jpg";
    specularMapPath = "images/soundbox.jpg";
    Cube drawing_sound_box = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 4.0f, 1.0f);

    //  drawing_door-------------------------------------------------------------
    diffuseMapPath = "images/door.PNG";
    specularMapPath = "images/door.PNG";
    Cube door = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    //  drawing_door-------------------------------------------------------------
    diffuseMapPath = "images/leaf.jpg";
    specularMapPath = "images/leaf.jpg";
    Cube leaf_texture = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    //  sunroof-------------------------------------------------------------
    diffuseMapPath = "images/sunroof.png";
    specularMapPath = "images/sunroof.png";
    Cube sunroof = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);


    //  sunroof-------------------------------------------------------------
    diffuseMapPath = "images/glass.jpg";
    specularMapPath = "images/glass.jpg";
    Cube glass = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    //  sunroof-------------------------------------------------------------
    diffuseMapPath = "images/green_wood.jpg";
    specularMapPath = "images/green_wood.jpg";
    Cube green = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    
    // drawing_cupboard-----------------------------------------------------------------
    diffuseMapPath = "images/soundbox.jpg";
    specularMapPath = "images/soundbox.jpg";
    Cube drawing_cupboard = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 4.0f, 2.0f);

    // drawing_mat --------------------------------------------------------------------
    diffuseMapPath = "images/still.jpg";
    specularMapPath = "images/still.jpg";
    Cube drawing_mat = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // bed_room1_almari --------------------------------------------------------------------
    diffuseMapPath = "images/almari.jpg";
    specularMapPath = "images/almari.jpg";
    Cube drawing_almari = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 4.0f);

    // bed_room1_almari --------------------------------------------------------------------
    diffuseMapPath = "images/lotion.jpg";
    specularMapPath = "images/lotion.jpg";
    Cube door_left = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 4.0f);
    
    // bed_room1_almari --------------------------------------------------------------------
    diffuseMapPath = "images/almari_door_right.png";
    specularMapPath = "images/almari_door_right.png";
    Cube door_right = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 4.0f);

    // bed_room1_window ---------------------------------------------------------------------
    diffuseMapPath = "images/bed_room1_window.jpeg";
    specularMapPath = "images/bed_room1_window.jpeg";
    Cube drawing_bed_room1_window = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // bed_room1_bedtable ---------------------------------------------------------------------
    diffuseMapPath = "images/bed_side_table.jpg";
    specularMapPath = "images/bed_side_table.jpg";
    Cube drawing_bed_room1_bedtable = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 2.0f);

    // bed_room1_bed_sheet ---------------------------------------------------------------------
    diffuseMapPath = "images/bed_sheet.jpg";
    specularMapPath = "images/bed_sheet.jpg";
    Cube bed_sheet = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 3.0f);
    
    // bed_room1_bed_texture ---------------------------------------------------------------------
    diffuseMapPath = "images/bed_texture.PNG";
    specularMapPath = "images/bed_texture.PNG";
    Cube bed_texture = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 3.0f, 1.0f);
    
    // bed_room1_blanket_texture ---------------------------------------------------------------------
    diffuseMapPath = "images/clock.png";
    specularMapPath = "images/clock.png";
    Cube clock = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // bed_room1_blanket_texture ---------------------------------------------------------------------
    diffuseMapPath = "images/blanket_texture.jpg";
    specularMapPath = "images/blanket_texture.jpg";
    Cube blanket_texture = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 2.0f, 2.0f);

    // bed_room1_bed_pillow ---------------------------------------------------------------------
    diffuseMapPath = "images/bed_pillow.jpg";
    specularMapPath = "images/bed_pillow.jpg";
    Cube bed_pillow = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 3.0f);

    // dining_frize ---------------------------------------------------------------------
    diffuseMapPath = "images/floor.jpg";
    specularMapPath = "images/floor.jpg";
    Cube floor = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    //unsigned int textureID = dining_frize.getTextureID();  // Assuming Cube has this method

     // dining_frize ---------------------------------------------------------------------
    diffuseMapPath = "images/still.jpg";
    specularMapPath = "images/still.jpg";
    Cube dining_frize = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);


    // dining_frize ---------------------------------------------------------------------
    diffuseMapPath = "images/dining_window.jpg";
    specularMapPath = "images/dining_window.jpg";
    Cube dining_window = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);


    // kitchen_surface ---------------------------------------------------------------------
    diffuseMapPath = "images/kitchen_surface.jpg";
    specularMapPath = "images/kitchen_surface.jpg";
    Cube kitchen_surface = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 8.0f, 3.0f);

    // kitchen_surface_top ---------------------------------------------------------------------
    diffuseMapPath = "images/kitchen_surface_top.jpg";
    specularMapPath = "images/kitchen_surface_top.jpg";
    Cube kitchen_surface_top = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 5.0f, 2.0f);
    
    // kitchen_stove ---------------------------------------------------------------------
    diffuseMapPath = "images/stove.PNG";
    specularMapPath = "images/stove.PNG";
    Cube kitchen_stove = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // kitchen_back ---------------------------------------------------------------------
    diffuseMapPath = "images/kitchen_back_texture.jpg";
    specularMapPath = "images/kitchen_back_texture.jpg";
    Cube kitchen_back_texture = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // kitchen_cupboard ---------------------------------------------------------------------
    diffuseMapPath = "images/kitchen_cupboard.PNG";
    specularMapPath = "images/kitchen_cupboard.PNG";
    /*diffuseMapPath = "images/kitchen_back_texture.jpg";
    specularMapPath = "images/kitchen_back_texture.jpg";*/
    Cube kitchen_cupboard1 = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 5.0f, 1.0f);
    Cube kitchen_cupboard2 = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 2.0f, 1.0f);

    // bed_room2_bookshelf
    diffuseMapPath = "images/book.jpg";
    specularMapPath = "images/book.jpg";
    Cube bookshelf = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // roof
    diffuseMapPath = "images/roof.jpg";
    specularMapPath = "images/roof.jpg";
    Cube roof = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 12.0f, 8.0f);

    // chair
    diffuseMapPath = "images/chair_leather.jpg";
    specularMapPath = "images/chair_leather.jpg";
    Cube chair_texture = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // chair
    diffuseMapPath = "images/divan.jpg";
    specularMapPath = "images/divan.jpg";
    Cube divan = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // chair
    diffuseMapPath = "images/chair_texture2.jpg";
    specularMapPath = "images/chair_texture2.jpg";
    Cube chair_texture2 = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // dining_table_texture
    diffuseMapPath = "images/dining_table_texture.PNG";
    specularMapPath = "images/dining_table_texture.PNG";
    Cube dining_table_texture = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 2.0f, 1.0f);

    // dressing_right
    diffuseMapPath = "images/dressing_right.PNG";
    specularMapPath = "images/dressing_right.PNG";
    Cube dressing_right = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // dressing_bottom
    diffuseMapPath = "images/dressing_bottom.PNG";
    specularMapPath = "images/dressing_bottom.PNG";
    Cube dressing_bottom = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // dressing_mirror
    diffuseMapPath = "images/mirror.jpg";
    specularMapPath = "images/mirror.jpg";
    Cube mirror = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // bathroom_top
    diffuseMapPath = "images/bathroom_top.PNG";
    specularMapPath = "images/bathroom_top.PNG";
    Cube bathroom_top = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // bathroom_top
    diffuseMapPath = "images/bathroom_door.PNG";
    specularMapPath = "images/bathroom_door.PNG";
    Cube bathroom_door = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    
    // bathroom_tiles
    diffuseMapPath = "images/bathroom_tiles.PNG";
    specularMapPath = "images/bathroom_tiles.PNG";
    Cube bathroom_tiles = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 4.0f, 2.0f);

    // bathroom_toilet
    diffuseMapPath = "images/toilet.PNG";
    specularMapPath = "images/toilet.PNG";
    Cube bathroom_toilet = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // bathroom_toilet
    diffuseMapPath = "images/sofa_back.jpg";
    specularMapPath = "images/sofa_back.jpg";
    Cube curve_sofa = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 2.0f, 1.0f);

    // bathroom_toilet
    diffuseMapPath = "images/sofa.jpg";
    specularMapPath = "images/sofa.jpg";
    Cube dd = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    glm::mat4 alTogether = glm::mat4(1.0f);  // Identity matrix as a starting point

    //cone texture
    string diffuseMapPath9 = "images/sponge.jpg";
    string specularMapPath9 = "images/sponge.jpg";
    unsigned int diffMap9 = loadTexture(diffuseMapPath9.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap9 = loadTexture(specularMapPath9.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    // You can modify the 'alTogether' matrix with translations, rotations, scalings, etc.
    // Example: Here we add translation and scale transformations
    alTogether = glm::translate(alTogether, glm::vec3(2.0f, -1.0f, 0.5f));  // Translation
    alTogether = glm::scale(alTogether, glm::vec3(0.05f, 0.05f, 0.05f));

    //cone texture
    string diffuseMapPath5 = "images/lipstick.jpg";
    string specularMapPath5 = "images/lipstick.jpg";
    unsigned int diffMap5 = loadTexture(diffuseMapPath5.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap5 = loadTexture(specularMapPath5.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);

    //sphere texture
    string diffuseMapPath8 = "images/table_top2.jpg";
    string specularMapPath8 = "images/table_top2.jpg";
    unsigned int diffMap8 = loadTexture(diffuseMapPath8.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap8 = loadTexture(specularMapPath8.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);

    //cylinder texture
    string diffuseMapPath7 = "images/sofa_back.jpg";
    string specularMapPath7 = "images/sofa_back.jpg";
    unsigned int diffMap7 = loadTexture(diffuseMapPath7.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap7 = loadTexture(specularMapPath7.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    //cylinder texture
    string diffuseMapPath6 = "images/chair_leather.jpg";
    string specularMapPath6 = "images/chair_leather .jpg";
    unsigned int diffMap6 = loadTexture(diffuseMapPath6.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap6 = loadTexture(specularMapPath6.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);


    // drawing_point_light
    Sphere clock_bell = Sphere();
    //Cylinder drawingMat(0.6f, 0.3f, 0.01f, 36, 1); // Radius = 0.3, height = 0.01, 36 sectors

    //ourShader.use();
    //lightingShader.use();

    

    
    

    Curve buck(Bucket);
    bucket = &buck;
    Curve buck2(Backrest);
    backrest = &buck2;
    Curve S(Seat);
    seat = &S;
    Curve buck3(Backrest2);
    backrest2 = &buck3;

    Cone cone = Cone();

    pointlight2.turnOff();

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
        // --------------------
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // be sure to activate shader when setting uniforms/drawing objects
        lightingShader.use();
        lightingShader.setVec3("viewPos", camera.Position);

        //// point light 1
        //pointlight1.setUpPointLight(lightingShader);
        //// point light 2
        //pointlight2.setUpPointLight(lightingShader);
        //// point light 3
        //pointlight3.setUpPointLight(lightingShader);
        //// point light 4
        //pointlight4.setUpPointLight(lightingShader);


        // *************************************DRAWING ROOM LIGHT********************************
        //drawing_light.setUpPointLight(lightingShader);

        
        // activate shader
        lightingShader.use();

        // pass projection matrix to shader (note that in this case it could change every frame)
        // glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        //glm::mat4 projection = glm::ortho(-2.0f, +2.0f, -1.5f, +1.5f, 0.1f, 100.0f);
        projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        lightingShader.setMat4("projection", projection);

        // camera/view transformation
        // glm::mat4 view = camera.GetViewMatrix();
        //glm::mat4 view = basic_camera.createViewMatrix();
        view = camera.GetViewMatrix();
        lightingShader.setMat4("view", view);

        // Modelling Transformation
        glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
        glm::mat4 trans4,trans3,translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;
        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X, translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X, scale_Y, scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
        lightingShader.setMat4("model", model);

        //glBindVertexArray(cubeVAO);
        //glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        //glDrawArrays(GL_TRIANGLES, 0, 36);

        /*bed(cubeVAO, lightingShader, model,0);
        bed(cubeVAO, lightingShader, model, 2);*/
        //drawingRoom(cubeVAO, lightingShader, model, lightingShaderWithTexture);
        
        /*building(cubeVAO, lightingShader, model);
        road(cubeVAO, lightingShader, model);*/

        

        // ********************************DRAWING ROOM************************************
        // Base transformation for lipstick
        glm::mat4 lipstickTransform = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, 0.0f));

        // Base transformation for the showpiece
        glm::mat4 baseTransform = glm::mat4(1.0f);
        baseTransform = glm::translate(baseTransform, glm::vec3(0.0f, 0.0f, 0.0f)); // Position the showpiece
        rotationAngle += deltaTime * rotationSpeed; // Increment angle over time
        drawShowpiece(lightingShaderWithTexture, baseTransform,diffMap6, specMap6, diffMap5, specMap5, rotationAngle);
        // Render the lipstick
        //drawLipstick(lightingShaderWithTexture, lipstickTransform, diffMap6, specMap6,diffMap5,specMap5);

        drawPedestal(lightingShaderWithTexture,diffMap7,specMap7);

        // Example call to draw the mirror with the stand
        drawMirrorWithStand(lightingShaderWithTexture, alTogether, floor, mirror, floor, 0.0f, 20.0f, -10.0f, 0.0f,6.0f);
        glm::mat4 als = glm::mat4(1.0f); // Identity matrix for transformation

        duvan(chair_texture, divan, lightingShaderWithTexture, als, diffMap6, specMap6);
        
        // Initial parameters for drawing
        float tx = 0.0f; // Initial x-position
        float ty = 20.0f; // Initial y-position
        float tz = 15.0f; // Initial z-position
        int depth = 3; // Maximum recursion depth
        float scaleFactor = 4.0f; // Scale factor for the cactus
        float branchAngle = 30.0f; // Angle of the branches

        // Function call to start drawing the cactus fractal
        drawFractalCactus(cubeVAO, lightingShaderWithTexture, alTogether, bathroom_tiles, leaf_texture, tx, ty, tz, depth, scaleFactor, branchAngle);
        f_waiting_floor(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_floor);
        f_waiting_wall(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_wall, door, sunroof,glass,green,ourShader);
        f_waiting_sofa(cubeVAO, lightingShader, model, lightingShaderWithTexture, sofa_top, sofa_foam, sofa_pillow);
        f_waiting_sofa1(cubeVAO, lightingShader, model, lightingShaderWithTexture, sofa_top, sofa_foam, sofa_pillow, 90.0, 0, 0, 0);
        f_waiting_table(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_table);
        //f_drawing_tv(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_tv, drawing_sound_box, drawing_cupboard);
        //f_waiting_mat(cylinderVAO, lightingShader, alTogether, lightingShaderWithTexture, drawingMat);
        //fan(cubeVAO, lightingShader, model, lightingShaderWithTexture, b, c, 0.0f, 0.0f, 0.0f);
        //f_drawing_window(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_window);
        
        //  **********************************BED ROOM 1*************************************
        f_room1_wall(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_wall, ourShader, clock_bell, clock, dd);
        f_room1_floor(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_floor);
        //f_bed_room1_bed(cubeVAO, lightingShader, model, lightingShaderWithTexture, bed_sheet, bed_pillow, bed_texture, blanket_texture);
        f_room1_almari(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_almari,door_right,door_left);
        
        //f_bed_room1_window(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_bed_room1_window);
        f_room1_bedtable(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_bed_room1_bedtable, ourShader);
        //fan(cubeVAO, lightingShader, model, lightingShaderWithTexture, b, c, 0.0f, 0.0f, -2.8f);
        
        // ***************************************DINING ROOM***************************************
        f_room2_wall(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_wall);
        f_room2_floor(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_floor);
        f_room2_table(cubeVAO, lightingShader, model, lightingShaderWithTexture, dining_table_texture,mirror);
        //f_dining_frize(cubeVAO, lightingShader, model, lightingShaderWithTexture, dining_frize);
        f_room2_window(cubeVAO, lightingShader, model, lightingShaderWithTexture, dining_window);
        //fan(cubeVAO, lightingShader, model, lightingShaderWithTexture, b, c, 3.0f, 0.0f, -2.5f);
        f_room2_light(cubeVAO, lightingShader, model, lightingShaderWithTexture, ourShader);
        chair(cubeVAO, lightingShaderWithTexture, model, chair_texture,dining_frize, 2.0f, 0.4f, -3.0f, chairRotation);
        chair(cubeVAO, lightingShaderWithTexture, model, chair_texture, dining_frize, 3.0f, 0.4f, -3.0f, chairRotation);
        chair(cubeVAO, lightingShaderWithTexture, model, chair_texture, dining_frize, 4.0f, 0.4f, -3.0f, chairRotation);
        

        chair2(cubeVAO, lightingShaderWithTexture, model, chair_texture2, dining_frize, -8.0f, 1.5f, 0.3f, -90.0);
        chair2(cubeVAO, lightingShaderWithTexture, model, chair_texture2, dining_frize, -11.0f, 1.5f, 0.3f, -90.0);
        chair2(cubeVAO, lightingShaderWithTexture, model, chair_texture2, dining_frize, -14.0f, 1.5f, 0.3f, -90.0);

        f_room2_table2(cubeVAO, lightingShader, model, lightingShaderWithTexture, dining_table_texture, mirror);
        // ****************************************KITCHEN******************************************
        //f_kitchen_surface(cubeVAO, lightingShader, model, lightingShaderWithTexture, kitchen_surface, kitchen_surface_top, kitchen_cupboard1, kitchen_cupboard2, kitchen_back_texture);
        //f_kitchen_stove(cubeVAO, lightingShader, model, lightingShaderWithTexture, kitchen_stove);
        
        // *****************************************BED ROOM 2**************************************
        f_room3_wall(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_wall);
        f_room3_floor(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_floor);
        //f_bed_room2_book(cubeVAO, lightingShader, model, lightingShaderWithTexture, bookshelf);
        //fan(cubeVAO, lightingShader, model, lightingShaderWithTexture, b, c, 3.5f, 0.0f, 0.5f);
        f_room3_light(cubeVAO, lightingShader, model, lightingShaderWithTexture, ourShader);
        //f_bed_room2_dressing_table(cubeVAO, lightingShader, model, lightingShaderWithTexture, dressing_right, dressing_bottom, dressing_mirror);
        //f_bed_room2_bed(cubeVAO, lightingShader, model, lightingShaderWithTexture, bed_sheet, bed_pillow, bed_texture, blanket_texture);
        
        // ********************************************ROOF******************************************
        f_roof(cubeVAO, lightingShader, model, lightingShaderWithTexture, roof);

        // ***********************************************BATHROOM***********************************
        //f_bathroom_wall(cubeVAO, lightingShader, model, lightingShaderWithTexture, bathroom_top, bathroom_door, bathroom_tiles);
        //f_bathroom_light(cubeVAO, lightingShader, model, lightingShaderWithTexture, ourShader);
        //f_toilet(cubeVAO, lightingShader, model, lightingShaderWithTexture, bathroom_toilet);
        f_bathroom_bucket(cubeVAO, lightingShader, model, lightingShaderWithTexture, curve_sofa);

        // also draw the lamp object(s)
        ourShader.use();
        ourShader.setMat4("projection", projection);
        ourShader.setMat4("view", view);

        // we now draw as many light bulbs as we have point lights.
        glBindVertexArray(lightCubeVAO);
        //for (unsigned int i = 0; i < 4; i++)
        //{
        //    model = glm::mat4(1.0f);
        //    model = glm::translate(model, pointLightPositions[i]);
        //    model = glm::scale(model, glm::vec3(0.2f)); // Make it a smaller cube
        //    ourShader.setMat4("model", model);
        //    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
        //    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        //    //glDrawArrays(GL_TRIANGLES, 0, 36);
        //}

        ////texture
        //glm::mat4 modelTexture = glm::mat4(1.0f);
        //glm::mat4 translate = glm::mat4(1.0f);
        //glm::mat4 scale = glm::mat4(0.5f);

        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X+1, translate_Y+1, translate_Z + 1));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X * 0.5, scale_Y * 0.5, scale_Z * 0.5));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
        
        lightingShaderWithTexture.use();
        lightingShaderWithTexture.setVec3("viewPos", camera.Position);
        lightingShaderWithTexture.setMat4("view", view);
        lightingShaderWithTexture.setMat4("projection", projection);

        lightingShaderWithTexture.use();
        // point light 1
        pointlight1.setUpPointLight(lightingShaderWithTexture);
        // point light 2
        pointlight2.setUpPointLight(lightingShaderWithTexture);
        // point light 3
        pointlight3.setUpPointLight(lightingShaderWithTexture);
        // point light 4
        pointlight4.setUpPointLight(lightingShaderWithTexture);

        //  *****************************DRAWING_LIGHT*****************************
        drawing_light.setUpPointLight(lightingShaderWithTexture);
        //  ******************************BED ROOM 1********************************
        bed_room1_light.setUpPointLight(lightingShaderWithTexture);
        //   ******************************DINING ROOM******************************
        dining_light.setUpPointLight(lightingShaderWithTexture);
        //   *******************************BED ROOM2*******************************
        bed_room2_light.setUpPointLight(lightingShaderWithTexture);
        //   ******************************BATHROOM LIGHT*********************************
        bathroom_light.setUpPointLight(lightingShaderWithTexture);

        /*diffuseMapPath = "images/emoji.png";
        specularMapPath = "images/emoji.png";
        diffMap = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
        specMap = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
        Cube cube = Cube(diffMap, specMap, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
        cube.drawCubeWithTexture(lightingShaderWithTexture, model);*/

        

        //cone
        translateMatrix = glm::translate(identityMatrix, glm::vec3(3.0f, 0.3f, -1.8f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(0.08f, 0.08f, 0.08f));
        trans3 = glm::translate(identityMatrix, glm::vec3(0.0f, 0.0f, bathroom_door_translate*7.0f));
        trans4 = glm::translate(identityMatrix, glm::vec3(0.0f, spongeup, 0.0f));
        model = translateMatrix * scaleMatrix*trans3*trans4;
        

        glm::vec3 color = glm::vec3(0.494f, 0.514f, 0.541f);
        Cone cone = Cone(0.4f, 1.0f, 36, color, color, color, 32.0f);
        cone.drawConeTexture(lightingShaderWithTexture, model, diffMap9, specMap9);
        //cone
        translateMatrix = glm::translate(identityMatrix, glm::vec3(2.9f, 0.3f, -1.8f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(0.08f, 0.08f, 0.08f));
        trans3 = glm::translate(identityMatrix, glm::vec3(0.0f, 0.0f, bathroom_door_translate * 7.0f));
        trans4 = glm::translate(identityMatrix, glm::vec3(0.0f, spongeup, 0.0f));
        model = translateMatrix * scaleMatrix * trans3 * trans4;


        glm::vec3 color2 = glm::vec3(0.494f, 0.514f, 0.541f);
        Cone cone2 = Cone(0.4f, 1.0f, 36, color2, color2, color2, 32.0f);
        cone2.drawConeTexture(lightingShaderWithTexture, model, diffMap9, specMap9);

        //sphere
        
        translateMatrix = glm::translate(identityMatrix, glm::vec3(-1.0f, 0.5f, -0.8f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(0.3f, 0.3f, 0.3f));
        model = translateMatrix * scaleMatrix;
        color = glm::vec3(0.494f, 0.514f, 0.541f);
        Sphere sphere = Sphere(0.4f, 1.0f, 36, color, color, color, 32.0f);
        sphere.drawSphereTexture(lightingShaderWithTexture, model, diffMap8, specMap8);

        // Chandelier position
        float chandelierX = 0.0f;
        float chandelierY = 20.0f;
        float chandelierZ = 0.0f;

        // Transformation matrix
        glm::mat4 alTogether = glm::mat4(1.0f);

        // Call the chandelier function
        chandelier(cubeVAO, lightingShader, alTogether, lightingShaderWithTexture, dining_frize, sphere, diffMap8, specMap8, 0.0f, 1.0f, 0.0f);

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    glDeleteVertexArrays(1, &cubeVAO);
    glDeleteVertexArrays(1, &lightCubeVAO);
    glDeleteBuffers(1, &cubeVBO);
    glDeleteBuffers(1, &cubeEBO);

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}

void shaderActivate(Shader& shader)
{
    shader.use();
    shader.setVec3("viewPos", camera.Position);
    shader.setMat4("view", view);
    shader.setMat4("projection", projection);
}


void drawCube(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 model = glm::mat4(1.0f), float r = 1.0f, float g = 1.0f, float b = 1.0f)
{
    lightingShader.use();

    lightingShader.setVec3("material.ambient", glm::vec3(r, g, b));
    lightingShader.setVec3("material.diffuse", glm::vec3(r, g, b));
    lightingShader.setVec3("material.specular", glm::vec3(0.5f, 0.5f, 0.5f));
    lightingShader.setFloat("material.shininess", 32.0f);

    lightingShader.setMat4("model", model);

    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
}

// *******************waiting ROOM*****************************

void drawLipstick(Shader& shader, glm::mat4 baseTransform, unsigned int diffuseMap, unsigned int specularMap, unsigned int diffuseMapcone, unsigned int specularMapcone) {
    glm::mat4 model;
    glm::mat4 translateMatrix, scaleMatrix;

    // Apply uniform scaling to the entire lipstick structure
    float uniformScale = 0.1f; // Scale the lipstick to half its size
    baseTransform = glm::scale(baseTransform, glm::vec3(uniformScale, uniformScale, uniformScale));
    // Apply uniform translation to the entire lipstick structure
    glm::vec3 uniformTranslation = glm::vec3(30.0f, 3.5f, -18.5f); // Translate to (1.0, 0.5, -2.0)
    baseTransform = glm::translate(baseTransform, uniformTranslation);
    // Create cylinder and cone objects
    Cylinder baseCylinder(0.3f, 0.5f, 36);    // Base cylinder (container)
    Cylinder innerCylinder(0.15f, 0.4f, 36); // Inner cylinder (lipstick body)
    Cone lipstickTip(0.15f, 0.2f, 36);       // Lipstick tip (cone shape)

    // Base Cylinder (Container)
    translateMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, -0.2f, 0.0f)); // Lower position
    scaleMatrix = glm::scale(glm::mat4(1.0f), glm::vec3(0.8f, 1.0f, 0.8f)); // No scaling
    model = baseTransform * translateMatrix * scaleMatrix;
    shader.setMat4("model", model);
    baseCylinder.drawCylinderTexture(shader, model, diffuseMap, specularMap);

    // Inner Cylinder (Lipstick Body)
    translateMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.25f, 0.0f)); // Raise it above the container
    model = baseTransform * translateMatrix;
    shader.setMat4("model", model);
    innerCylinder.drawCylinderTexture(shader, model, diffuseMap, specularMap);

    // Lipstick Tip (Cone)
    translateMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.45f, 0.0f)); // Position on top of the inner cylinder
    model = baseTransform * translateMatrix;
    shader.setMat4("model", model);
    lipstickTip.drawConeTexture(shader, model, diffuseMapcone, specularMapcone);
}

void drawPedestal(Shader& lightingShaderWithTexture, unsigned int diffuseMap, unsigned int specularMap) {
    // Step dimensions
    float stepHeight = 0.02f;
    float stepRadius[] = { 0.3f, 0.27f, 0.25f };
    glm::vec3 stepPositions[] = {
        glm::vec3(3.0f, 0.02f, 0.5f), // Bottom step
        glm::vec3(3.0f, 0.04f, 0.5f), // Middle step
        glm::vec3(3.0f,  0.06f, 0.5f)  // Top step
    };

    // Create Cylinder objects for each step
    Cylinder step1(stepRadius[0], stepHeight, 36);
    Cylinder step2(stepRadius[1], stepHeight, 36);
    Cylinder step3(stepRadius[2], stepHeight, 36);

    // Render each step with appropriate transformations
    glm::mat4 model;

    // Bottom step
    model = glm::mat4(1.0f);
    model = glm::translate(model, stepPositions[0]);
    step1.drawCylinderTexture(lightingShaderWithTexture, model, diffuseMap, specularMap);

    // Middle step
    model = glm::mat4(1.0f);
    model = glm::translate(model, stepPositions[1]);
    step2.drawCylinderTexture(lightingShaderWithTexture, model, diffuseMap, specularMap);

    // Top step
    model = glm::mat4(1.0f);
    model = glm::translate(model, stepPositions[2]);
    step3.drawCylinderTexture(lightingShaderWithTexture, model, diffuseMap, specularMap);
}


void duvan(Cube& cube2, Cube& cube3, Shader& lightingShaderWithTexture, glm::mat4 als, unsigned int diffuseMap, unsigned int specularMap) {
    glm::mat4 identityMatrix = glm::mat4(1.0f); // Initialize to identity matrix
    glm::mat4 translateMatrix, rotateYMatrix, scaleMatrix, model;

    // Apply global translation and scaling to the entire object
    glm::vec3 globalTranslation = glm::vec3(3.7f, 0.0f, -0.1f); // Move the object to (2.0, 0.0, -5.0)
    float uniformScale = 0.5f; // Scale the entire object to half its size
    als = glm::translate(als, globalTranslation); // Apply global translation
    als = glm::scale(als, glm::vec3(uniformScale, uniformScale, uniformScale)); // Apply global scaling

    // Create Cylinder Objects
    Cylinder pillowCylinder(0.1f, 0.5f, 36);     // Pillow cylinder with radius 0.1 and height 0.5
    Cylinder supportCylinder(0.075f, 0.5f, 36);  // Back support cylinder with smaller radius

    // Base
    translateMatrix = glm::translate(identityMatrix, glm::vec3(0.00, 0.5, 0.0));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0f, .125, 1.7f));
    model = als * translateMatrix * scaleMatrix;
    lightingShaderWithTexture.setMat4("model", model);
    cube2.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Leg 1
    translateMatrix = glm::translate(identityMatrix, glm::vec3(0.00, 0.5, 0.0));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0f, -.35, 0.05f));
    model = als * translateMatrix * scaleMatrix;
    lightingShaderWithTexture.setMat4("model", model);
    cube2.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Leg 2
    translateMatrix = glm::translate(identityMatrix, glm::vec3(0.00, 0.5, 1.7 - 0.05));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0f, -.35, 0.05f));
    model = als * translateMatrix * scaleMatrix;
    lightingShaderWithTexture.setMat4("model", model);
    cube2.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Seat
    translateMatrix = glm::translate(identityMatrix, glm::vec3(0.00, 0.75, -0.1));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0f, -.125, 1.85f));
    model = als * translateMatrix * scaleMatrix;
    lightingShaderWithTexture.setMat4("model", model);
    cube3.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Back
    translateMatrix = glm::translate(identityMatrix, glm::vec3(0.00, 1.05, -0.45));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(45.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0f, -.05, 0.7f));
    model = als * translateMatrix * rotateYMatrix * scaleMatrix;
    lightingShaderWithTexture.setMat4("model", model);
    cube2.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Back Seat
    translateMatrix = glm::translate(identityMatrix, glm::vec3(0.00, 1.05, -0.40));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(45.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0f, -.05, 0.6f));
    model = als * translateMatrix * rotateYMatrix * scaleMatrix;
    lightingShaderWithTexture.setMat4("model", model);
    cube3.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Front
    translateMatrix = glm::translate(identityMatrix, glm::vec3(0.00, 1.05, 2.20));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(135.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0f, -.05, 0.7f));
    model = als * translateMatrix * rotateYMatrix * scaleMatrix;
    lightingShaderWithTexture.setMat4("model", model);
    cube2.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Front Seat
    translateMatrix = glm::translate(identityMatrix, glm::vec3(0.00, 1.05, 2.15));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(135.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.0f, -.05, 0.6f));
    model = als * translateMatrix * rotateYMatrix * scaleMatrix;
    lightingShaderWithTexture.setMat4("model", model);
    cube3.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Pillow 1
    translateMatrix = glm::translate(identityMatrix, glm::vec3(0.5, 0.85, 0.0));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(90.f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.1, 1.5, 1.1));
    model = als * translateMatrix * rotateYMatrix * scaleMatrix;
    lightingShaderWithTexture.setMat4("model", model);
    pillowCylinder.drawCylinderTexture(lightingShaderWithTexture, model, diffuseMap, specularMap);

    // Pillow 2
    translateMatrix = glm::translate(identityMatrix, glm::vec3(0.5, 0.85, 1.75));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(90.f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(1.1, 1.5, 1.1));
    model = als * translateMatrix * rotateYMatrix * scaleMatrix;
    lightingShaderWithTexture.setMat4("model", model);
    pillowCylinder.drawCylinderTexture(lightingShaderWithTexture, model, diffuseMap, specularMap);

    // Back Support 1
    translateMatrix = glm::translate(identityMatrix, glm::vec3(0.5, 1.05, -0.45));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(90.f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(0.175, 1.5, 0.175));
    model = als * translateMatrix * rotateYMatrix * scaleMatrix;
    lightingShaderWithTexture.setMat4("model", model);
    supportCylinder.drawCylinderTexture(lightingShaderWithTexture, model, diffuseMap, specularMap);

    // Back Support 2
    translateMatrix = glm::translate(identityMatrix, glm::vec3(0.5, 1.05, 2.15));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(90.f), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(0.175, 1.5, 0.175));
    model = als * translateMatrix * rotateYMatrix * scaleMatrix;
    lightingShaderWithTexture.setMat4("model", model);
    supportCylinder.drawCylinderTexture(lightingShaderWithTexture, model, diffuseMap, specularMap);
}



void chandelier(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube, Sphere& sphere, unsigned int diffMap8, unsigned int specMap8, float x, float y, float z)
{
    // Central Hub (Reduced Size)
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 rotation = glm::mat4(1.0f);

    // Apply rotation to the entire chandelier
    rotation = glm::rotate(glm::mat4(1.0f), glm::radians(rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    scale = glm::scale(model, glm::vec3(0.02f, 0.02f, 0.02f)); // Smaller central hub
    translate = glm::translate(model, glm::vec3(-0.1f, 0.0f, -0.1f));
    translate3 = glm::translate(model, glm::vec3(x, y, z));
    model = alTogether * translate3 * rotation * scale * translate;
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Arms of the chandelier (Shortened and Slimmer)
    float armLength = 0.2f; // Shorter arms
    float armWidth = 0.01f; // Thinner arms
    float armHeight = 0.01f; // Slimmer arms
    int armCount = 6; // Number of arms
    for (int i = 0; i < armCount; ++i)
    {
        float angle = glm::radians(360.0f / armCount * i);
        glm::mat4 armRotation = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(0.0f, 1.0f, 0.0f));
        scale = glm::scale(glm::mat4(1.0f), glm::vec3(armLength, armHeight, armWidth));
        translate = glm::translate(glm::mat4(1.0f), glm::vec3(armLength / 2 - 0.1f, 0.0f, 0.0f)); // Move arm outward
        model = alTogether * translate3 * rotation * armRotation * translate * scale;
        cube.drawCubeWithTexture(lightingShaderWithTexture, model);

        // Bulbs at the end of the arms (Smaller)
        glm::mat4 bulbTranslate = glm::translate(glm::mat4(1.0f), glm::vec3(armLength, 0.0f, 0.0f));
        glm::mat4 bulbScale = glm::scale(glm::mat4(1.0f), glm::vec3(0.1f, 0.1f, 0.1f)); // Reduced bulb size
        model = alTogether * translate3 * rotation * armRotation * bulbTranslate * bulbScale;
        sphere.drawSphereTexture(lightingShaderWithTexture, model, diffMap8, specMap8);
    }

    // Vertical pole connecting to the ceiling (Slimmer and Shorter)
    glm::mat4 poleScale = glm::scale(glm::mat4(1.0f), glm::vec3(0.015f, 0.5f, 0.015f)); // Slimmer pole
    glm::mat4 poleTranslate = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, 0.0f)); // Shortened pole height
    model = alTogether * translate3 * rotation * poleTranslate * poleScale;
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
}


void drawFractalCactus(unsigned int& cubeVAO, Shader& lightingShaderWithTexture, glm::mat4 alTogether,
    Cube& base_texture, Cube& cactus_texture,
    float tx, float ty, float tz, int depth, float scaleFactor, float branchAngle) {
    if (depth == 0) return; // Stop recursion at depth 0

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 temp = glm::mat4(1.0f);

    // Draw base (pot as a cube)
    if (depth ==3) { // Base case for the ground or pot
        temp = glm::translate(glm::mat4(1.0f), glm::vec3(tx-0.06, ty, tz-0.06)) *
            glm::scale(glm::mat4(1.0f), glm::vec3(1.0f * scaleFactor, 0.5f * scaleFactor, 1.0f * scaleFactor));
        model = alTogether * temp;
        base_texture.drawCubeWithTexture(lightingShaderWithTexture, model);
    }

    // Draw the main cactus trunk segment
    temp = glm::translate(glm::mat4(1.0f), glm::vec3(tx, ty + 0.5f * scaleFactor, tz)) *
        glm::scale(glm::mat4(1.0f), glm::vec3(0.2f * scaleFactor, 1.0f * scaleFactor, 0.2f * scaleFactor));
    model = alTogether * temp;
    cactus_texture.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Branching logic for arms
    float newTx = tx;
    float newTy = ty + 0.8f * scaleFactor;
    float newTz = tz;

    // Left arm
    if (depth % 2 == 0) { // Alternate arm generation
        glm::mat4 rotateLeft = glm::rotate(glm::mat4(1.0f), glm::radians(branchAngle), glm::vec3(0.0f, 0.0f, 1.0f));
        glm::mat4 scaleLeft = glm::scale(glm::mat4(1.0f), glm::vec3(0.2f * scaleFactor, 0.6f * scaleFactor, 0.2f * scaleFactor));
        temp = glm::translate(glm::mat4(1.0f), glm::vec3(newTx - 0.1f * scaleFactor, newTy, newTz)) * rotateLeft * scaleLeft;
        model = alTogether * temp;
        cactus_texture.drawCubeWithTexture(lightingShaderWithTexture, model);
        drawFractalCactus(cubeVAO, lightingShaderWithTexture, alTogether, base_texture, cactus_texture,
            newTx - 0.3f * scaleFactor, newTy + 0.3f * scaleFactor, newTz, depth - 1, scaleFactor * 0.7f, branchAngle);
    }

    // Right arm
    if (depth % 2 == 1) {
        glm::mat4 rotateRight = glm::rotate(glm::mat4(1.0f), glm::radians(-branchAngle), glm::vec3(0.0f, 0.0f, 1.0f));
        glm::mat4 scaleRight = glm::scale(glm::mat4(1.0f), glm::vec3(0.2f * scaleFactor, 0.6f * scaleFactor, 0.2f * scaleFactor));
        temp = glm::translate(glm::mat4(1.0f), glm::vec3(newTx + 0.1f * scaleFactor, newTy, newTz)) * rotateRight * scaleRight;
        model = alTogether * temp;
        cactus_texture.drawCubeWithTexture(lightingShaderWithTexture, model);
        drawFractalCactus(cubeVAO, lightingShaderWithTexture, alTogether, base_texture, cactus_texture,
            newTx + 0.3f * scaleFactor, newTy + 0.3f * scaleFactor, newTz, depth - 1, scaleFactor * 0.7f, branchAngle);
    }

    // Continue the main trunk upwards
    drawFractalCactus(cubeVAO, lightingShaderWithTexture, alTogether, base_texture, cactus_texture,
        tx, ty + 1.0f * scaleFactor, tz, depth - 1, scaleFactor * 0.9f, branchAngle);
}



void f_waiting_sofa(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& sofa_top, Cube& sofa_foam, Cube& sofa_pillow)
{
    float baseHeight = 0.2f;
    float width = 1.2f;
    float length = 0.4f;
    
    // base
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0.0, -0.5));
    translate2 = glm::translate(model, glm::vec3(width-0.5- 0.5, 0.0, -1.0));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.546, 0.335, 0.316);
    shaderActivate(lightingShaderWithTexture);
    sofa_top.drawCubeWithTexture(lightingShaderWithTexture, model);
    
    // foam
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.3, 0.03, 0.3));
    translate2 = glm::translate(model, glm::vec3(width - 0.19 - 0.5, 0.2, -1.0+0.05));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.3, 0.03, 0.3));
    translate2 = glm::translate(model, glm::vec3(width - 0.195-0.3 - 0.5, 0.2, -1.0 + 0.05));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.3, 0.03, 0.3));
    translate2 = glm::translate(model, glm::vec3(width - 0.2-0.3-0.3 - 0.5, 0.2, -1.0 + 0.05));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);

    // back
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, 0.1));
    translate2 = glm::translate(model, glm::vec3(width - 0.5 - 0.5, 0.2, -1.15));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    // right
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.12, baseHeight, 0.4));
    translate2 = glm::translate(model, glm::vec3(width+0.04 - 0.5, 0.2, -1.0));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    // left
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.12, baseHeight, 0.4));
    translate2 = glm::translate(model, glm::vec3(width - 1.04 - 0.5, 0.2, -1.0));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    // left pillo
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 0.1, 0.01));
    translate2 = glm::translate(model, glm::vec3(width-0.8 - 0.5, 0.23, -1.09));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);
    sofa_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);

    // right pillo
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 0.1, 0.01));
    translate2 = glm::translate(model, glm::vec3(width - 0.19 - 0.5, 0.23, -1.09));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);
    sofa_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

}


void f_waiting_sofa1(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& sofa_top, Cube& sofa_foam, Cube& sofa_pillow, float rotateY, float tx, float ty, float tz)
{
    float baseHeight = 0.2f;
    float width = 1.2f;
    float length = 0.4f;

    // base
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0.0, -0.5));
    translate2 = glm::translate(model, glm::vec3(width - 0.5 - 0.5, 0.0, -1.0));
    glm::mat4 rotateYMatrix = glm::rotate(model, glm::radians(rotateY), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.546, 0.335, 0.316);
    shaderActivate(lightingShaderWithTexture);
    sofa_top.drawCubeWithTexture(lightingShaderWithTexture, model);

    // foam
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.3, 0.03, 0.3));
    translate2 = glm::translate(model, glm::vec3(width - 0.19 - 0.5, 0.2, -1.0 + 0.05));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.3, 0.03, 0.3));
    translate2 = glm::translate(model, glm::vec3(width - 0.195 - 0.3 - 0.5, 0.2, -1.0 + 0.05));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.3, 0.03, 0.3));
    translate2 = glm::translate(model, glm::vec3(width - 0.2 - 0.3 - 0.3 - 0.5, 0.2, -1.0 + 0.05));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);

    // back
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, 0.1));
    translate2 = glm::translate(model, glm::vec3(width - 0.5 - 0.5, 0.2, -1.15));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    // right
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.12, baseHeight, 0.4));
    translate2 = glm::translate(model, glm::vec3(width + 0.04 - 0.5, 0.2, -1.0));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    // left
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.12, baseHeight, 0.4));
    translate2 = glm::translate(model, glm::vec3(width - 1.04 - 0.5, 0.2, -1.0));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    // left pillo
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 0.1, 0.01));
    translate2 = glm::translate(model, glm::vec3(width - 0.8 - 0.5, 0.23, -1.09));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);
    sofa_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);

    // right pillo
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 0.1, 0.01));
    translate2 = glm::translate(model, glm::vec3(width - 0.19 - 0.5, 0.23, -1.09));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);
    sofa_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

}


void f_waiting_table(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_table)
{
    float baseHeight = 0.02;
    float width = 0.7;
    float length = 0.4;

    // base
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0.0, -0.5));
    translate2 = glm::translate(model, glm::vec3(width - 0.5, 0.2, -0.2));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.1, 0.1, 0.1);
    shaderActivate(lightingShaderWithTexture);
    drawing_table.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

    // leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.2, 0.05));
    translate2 = glm::translate(model, glm::vec3(width-0.175, 0.0, -0.35));
    model = alTogether * translate2 * scale * translate;
    drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);

    // leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.2, 0.05));
    translate2 = glm::translate(model, glm::vec3(width - 0.175, 0.0, -0.025));
    model = alTogether * translate2 * scale * translate;
    drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);

    // leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.2, 0.05));
    translate2 = glm::translate(model, glm::vec3(width - 0.825, 0.0, -0.35));
    model = alTogether * translate2 * scale * translate;
    drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);

    // leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.2, 0.05));
    translate2 = glm::translate(model, glm::vec3(width - 0.825, 0.0, -0.025));
    model = alTogether * translate2 * scale * translate;
    drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);
    
}


void f_waiting_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube &cube)
{
    float baseHeight = 0.01;
    float width = 3.0;
    float length = 3.0;

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width+1, baseHeight, length+1));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}

/*
void f_waiting_mat(unsigned int& cylinderVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cylinder& drawing_mat)
{
    float baseHeight = 0.01f;  // Thickness of the mat
    float radius = 0.6f;       // Radius of the circular mat

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    // Scale the cylinder to match the mat's size
    scale = glm::scale(model, glm::vec3(radius, baseHeight, radius));

    // Position the mat above the table
    translate = glm::translate(model, glm::vec3(0.0f, 1.01f, 0.5f)); // Slightly above the table to avoid z-fighting

    // Combine transformations
    model = alTogether * translate * scale;

    // Render the circular mat with texture
    shaderActivate(lightingShaderWithTexture);
    drawing_mat.drawCylinder(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}
*/


void f_waiting_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube, Cube& door,Cube& sunroof, Cube& glass, Cube& green, Shader& ourShader)
{
    float baseHeight = 1.5;
    float width = 0.01;
    float length = 3.0;
   
    // right
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width+0.1, baseHeight, length+0.1));
    translate = glm::translate(model, glm::vec3(1.5, 0.0, -1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);

    // left
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width+0.1, baseHeight, length+0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // back
    width = 3.0;
    length = 0.01;
    //model = glm::mat4(1.0f);
    //scale = glm::scale(model, glm::vec3(width, baseHeight, length + 0.1));
    //translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    //model = alTogether * translate * scale;
    ////drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //

    // front top
    float baseHeight1 = baseHeight / 3;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight1, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 2*baseHeight1, 1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // front bottom
    float baseHeight2 = baseHeight1 * 2;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width-0.5, baseHeight2, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.0, 0.0, 1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    glass.drawCubeWithTexture(lightingShaderWithTexture, model);
    

    // back top
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight1, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 2 * baseHeight1, -1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // back bottom
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width - 0.5, baseHeight2, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // front door
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.5, baseHeight2, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 0.0, 1.5));
    glm::mat4 rotateYMatrix = glm::rotate(model, glm::radians(f_door), glm::vec3(0.0f, 1.0f, 0.0f));
    glm::mat4 translate_origin = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    model = alTogether * translate * rotateYMatrix * translate_origin * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    door.drawCubeWithTexture(lightingShaderWithTexture, model);
    /*    // back door
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.5, baseHeight2, length + 0.1));
    translate = glm::translate(model, glm::vec3(width-2.0, 0.0, -1.5));
    rotateYMatrix = glm::rotate(model, glm::radians(b1_door), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate * rotateYMatrix * translate_origin * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    door.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(ourShader);
    */

    // white sunroof above the door
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(3.0, 0.02, 0.4)); // Adjust size as needed
    translate = glm::translate(model, glm::vec3(-1.5, baseHeight2 + 0.05, 1.5)); // Position it just above the door
    glm::mat4 rotateXMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(20.0f), glm::vec3(1.0f, 0.0f, 0.0f)); // Rotate 30 degrees downward
    model = alTogether * translate * rotateXMatrix * scale;
    lightingShaderWithTexture.setVec3("color", glm::vec3(1.0f, 1.0f, 1.0f)); // Set color to white
    green.drawCubeWithTexture(lightingShaderWithTexture, model);

    // point light
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.15, 0.15));
    translate = glm::translate(model, glm::vec3(1.45f, 1.3f, 0.1f));
    model = alTogether * translate * scale;
    ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    shaderActivate(lightingShader);

    // Manually Placed Grid on Front Bottom Glass
    float gridThickness = 0.05; // Thickness of the grid bars
    float gridHeight = baseHeight2-0.1; // Height of the grid (same as front bottom glass)
    float gridWidth = length + 0.1; // Width of the grid (same as front bottom glass)

    // Vertical bars
    // Left bar
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(gridThickness, gridHeight, 0.02f));
    translate = glm::translate(model, glm::vec3(-1.0, 0.0, 1.6)); // Slightly in front of glass
    model = alTogether * translate * scale;
    green.drawCubeWithTexture(lightingShaderWithTexture, model);
    
    // Center bar
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(gridThickness, gridHeight, 0.02f));
    translate = glm::translate(model, glm::vec3(-0.5, 0.0, 1.6)); // Slightly in front of glass
    model = alTogether * translate * scale;
    green.drawCubeWithTexture(lightingShaderWithTexture, model);
    
    // Right bar
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(gridThickness, gridHeight, 0.02f));
    translate = glm::translate(model, glm::vec3(0.0, 0.0 , 1.6)); // Slightly in front of glass
    model = alTogether * translate * scale;
    green.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Center bar
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(gridThickness, gridHeight, 0.02f));
    translate = glm::translate(model, glm::vec3(0.5, 0.0, 1.6)); // Slightly in front of glass
    model = alTogether * translate * scale;
    green.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Right bar
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(gridThickness, gridHeight, 0.02f));
    translate = glm::translate(model, glm::vec3(1.0, 0.0, 1.6)); // Slightly in front of glass
    model = alTogether * translate * scale;
    green.drawCubeWithTexture(lightingShaderWithTexture, model);
    // Right bar
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(gridThickness, gridHeight, 0.02f));
    translate = glm::translate(model, glm::vec3(1.5, 0.0, 1.6)); // Slightly in front of glass
    model = alTogether * translate * scale;
    green.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Horizontal bars
    // Bottom bar
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.5, 0.5, 0.02f));
    translate = glm::translate(model, glm::vec3(-1.0, 0.0, 1.6)); // Slightly in front of glass
    model = alTogether * translate * scale;
    green.drawCubeWithTexture(lightingShaderWithTexture, model);

}


void f_waiting_tv(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_tv, Cube& drawing_sound_box, Cube& drawing_cupboard)
{
    float baseHeight = 0.4;
    float width = 0.6;
    float length = 0.02;

    // tv
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(0.0, 0.5, 1.475));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.0, 0.0, 0.0);
    shaderActivate(lightingShaderWithTexture);
    drawing_tv.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);

    // sound box
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, 0.05, length));
    translate = glm::translate(model, glm::vec3(0.0, 0.4, 1.475));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.0, 0.0, 0.0);
    drawing_sound_box.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);

    // cupboard
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, 0.2, length+0.15));
    translate = glm::translate(model, glm::vec3(0.0, 0.0, 1.475-0.15));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.0, 0.0, 0.0);
    drawing_cupboard.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);
}


void fan(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& b, Cube& c, float x, float y, float z)
{
    float bladel = 1.5;
    float bladew = 0.2;
    float bladeh = 0.01;

    // Center
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 scale2 = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.27, 0.3, 0.27));
    scale2 = glm::scale(model, glm::vec3(0.5, 0.5, 0.5));
    translate = glm::translate(model, glm::vec3(-0.67, 0.0, -0.4));
    translate2 = glm::translate(model, glm::vec3(0.0, 1.35, 0.0));
    translate3 = glm::translate(model, glm::vec3(x, y, z));
    model = alTogether * translate3 * translate2 * scale2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    shaderActivate(lightingShaderWithTexture);
    c.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bladel, bladeh, bladew));
    translate = glm::translate(model, glm::vec3(0.01, 0.0, 0.0));
    glm::mat4 rotateM = glm::rotate(model, glm::radians(45.0f + rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate3 * translate2 * scale2 * rotateM * scale * translate;
    shaderActivate(lightingShaderWithTexture);
    b.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bladel, bladeh, bladew));
    translate = glm::translate(model, glm::vec3(0, 0.01, 0.0));
    rotateM = glm::rotate(model, glm::radians(165.0f + rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate3 * translate2 * scale2 * rotateM * scale * translate;
    b.drawCubeWithTexture(lightingShaderWithTexture, model);


    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bladel, bladeh, bladew));
    translate = glm::translate(model, glm::vec3(0.01, 0.01, 0.0));
    rotateM = glm::rotate(model, glm::radians(285.0f + rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate3 * translate2 * scale2 * rotateM * scale * translate;
    b.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void f_waiting_window(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_window)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1, 0.7, 0.05));
    glm::mat4 rotateY = glm::rotate(model, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    translate = glm::translate(model, glm::vec3(-0.3, 0.5, -1.4));
    model = alTogether * rotateY * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    shaderActivate(lightingShaderWithTexture);
    drawing_window.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


// *************************ROOM 1*************************************
void f_room1_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube, Shader& ourShader, Sphere& clock_bell, Cube& clock, Cube& x)
{
    float baseHeight = 1.5;
    float width = 0.01;
    float length = 3.0;

    // right
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 rotateYMatrix = glm::rotate(model, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    scale = glm::scale(model, glm::vec3(width + 0.1, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(1.5, 0.0, -1.5));
    translate2 = glm::translate(model, glm::vec3(0.0, 0.0, -length));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);

    // left
    //model = glm::mat4(1.0f);
    //scale = glm::scale(model, glm::vec3(width + 0.1, baseHeight, length + 0.1));
    //translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    //model = alTogether * translate2 * rotateYMatrix * translate * scale;
    ////drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // back
    width = 3.3;
    length = 0.01;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // front top
    baseHeight = baseHeight / 3;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 2 * baseHeight, 1.5));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // front bottom
    baseHeight = baseHeight * 2;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width - 0.5, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.0, 0.0, 1.5));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // point light
    shaderActivate(ourShader);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.15, 0.15));
    translate = glm::translate(model, glm::vec3(1.0f, 1.3f, -3.1f));
    model = alTogether * translate * scale;
    ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    //shaderActivate(lightingShader);





    

}


void f_room1_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube)
{
    float baseHeight = 0.01;
    float width = 3.0;
    float length = 3.0;

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width + 1, baseHeight, length + 1));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -1.5));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}



void f_room1_almari(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& draw_almari,Cube& door_right, Cube& door_left)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    // Sides of the wardrobe
    glm::mat4 leftSide = glm::translate(model, glm::vec3(3.6, 0.0, 1.3)); // Position left wall
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.03, 1.0, 0.3));       // Thin vertical wall
    leftSide = alTogether * leftSide * scale;
    draw_almari.drawCubeWithTexture(lightingShaderWithTexture, leftSide);

    glm::mat4 rightSide = glm::translate(model, glm::vec3(3.0, 0.0, 1.3)); // Position right wall
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.03, 1.0, 0.3));        // Thin vertical wall
    rightSide = alTogether * rightSide * scale;
    draw_almari.drawCubeWithTexture(lightingShaderWithTexture, rightSide);

    // Top of the wardrobe
    glm::mat4 top = glm::translate(model, glm::vec3(3.0, 1.0, 1.3)); // Position top wall
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.6, 0.03, 0.3));  // Thin horizontal wall
    top = alTogether * top * scale;
    draw_almari.drawCubeWithTexture(lightingShaderWithTexture, top);

    // Bottom of the wardrobe
    glm::mat4 bottom = glm::translate(model, glm::vec3(3.0,0.0, 1.3)); // Position bottom wall
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.6, 0.03, 0.3));       // Thin horizontal wall
    bottom = alTogether * bottom * scale;
    draw_almari.drawCubeWithTexture(lightingShaderWithTexture, bottom);
    // mid1 of the wardrobe
    glm::mat4 mid1 = glm::translate(model, glm::vec3(3.0, 0.4, 1.3)); // Position bottom wall
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.6, 0.03, 0.3));       // Thin horizontal wall
    mid1 = alTogether * mid1 * scale;
    draw_almari.drawCubeWithTexture(lightingShaderWithTexture, mid1);
    // bot1 of the wardrobe
    glm::mat4 bot1 = glm::translate(model, glm::vec3(3.2, 0.4, 1.4)); // Position bottom wall
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.1, 0.2, 0.05));       // Thin horizontal wall
    bot1 = alTogether * bot1 * scale;
    door_left.drawCubeWithTexture(lightingShaderWithTexture, bot1);
    // bot2 of the wardrobe
    glm::mat4 bot2 = glm::translate(model, glm::vec3(3.4, 0.4, 1.4)); // Position bottom wall
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.1, 0.2, 0.05));       // Thin horizontal wall
    bot2 = alTogether * bot2 * scale;
    door_left.drawCubeWithTexture(lightingShaderWithTexture, bot2);

    // mid2 of the wardrobe
    glm::mat4 mid2 = glm::translate(model, glm::vec3(3.0, 0.7, 1.3)); // Position bottom wall
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.6, 0.03, 0.3));       // Thin horizontal wall
    mid2 = alTogether * mid2 * scale;
    draw_almari.drawCubeWithTexture(lightingShaderWithTexture, mid2);

    // Back of the wardrobe
    glm::mat4 back = glm::translate(model, glm::vec3(3.0, 0.0, 1.45)); // Position back wall
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.6, 1.0, 0.03));    // Thin back wall
    back = alTogether * back * scale;
    draw_almari.drawCubeWithTexture(lightingShaderWithTexture, back);

    // Front of the wardrobe
    glm::mat4 front = glm::translate(model, glm::vec3(3.0, 0.0, 1.3)); // Position back wall
    glm::mat4 rotateYMatrix = glm::rotate(model, glm::radians(b1_door), glm::vec3(0.0f, 1.0f, 0.0f));
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.6, 1.0, 0.03));    // Thin back wall
    front = alTogether * front * rotateYMatrix;
    front = front * scale;
    door_right.drawCubeWithTexture(lightingShaderWithTexture, front);

    

    shaderActivate(lightingShader);
}



void f_room1_bed(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& bed_sheet, Cube& bed_pillow, Cube& bed_texture, Cube& blanket_texture)
{
    float baseHeight = 0.3;
    float width = 1;
    float length = 2;
    float pillowWidth = 0.3;
    float pillowLength = 0.2;
    float blanketWidth = 0.8;
    float blanketLength = 0.7;
    float headHeight = 0.6;

    //base
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    translate3 = glm::translate(model, glm::vec3(-0.75, 0, -2.7));
    model = alTogether * translate3 * scale * translate;
    shaderActivate(lightingShaderWithTexture);
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    bed_texture.drawCubeWithTexture(lightingShaderWithTexture, model);

    //foam
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, baseHeight, 0));
    scale = glm::scale(model, glm::vec3(width, 0.06, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.804, 0.361, 0.361);
    bed_sheet.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

    //pillow 1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3((width / 2) - (0.1 + pillowWidth / 2), baseHeight + 1 * 0.06, (length / 2) - (0.025 + pillowWidth / 2)));
    scale = glm::scale(model, glm::vec3(pillowWidth, 0.04, pillowLength));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 1, 0.647, 0);
    shaderActivate(lightingShaderWithTexture);
    bed_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);

    //pillow 2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3((-width / 2) + (0.1 + pillowWidth / 2), baseHeight + 1 * 0.06, (length / 2) - (0.025 + pillowWidth / 2)));
    scale = glm::scale(model, glm::vec3(pillowWidth, 0.04, pillowLength));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 1, 0.647, 0);
    bed_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

    //blanket
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, baseHeight + 1 * 0.06, -(length / 2 - 0.025) + blanketLength / 2));
    scale = glm::scale(model, glm::vec3(blanketWidth, 0.015, blanketLength));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.541, 0.169, 0.886);
    blanket_texture.drawCubeWithTexture(lightingShaderWithTexture, model);

    //head
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, 0, (length / 2 - 0.02 / 2) + 0.02));
    scale = glm::scale(model, glm::vec3(width, headHeight, 0.06));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    bed_texture.drawCubeWithTexture(lightingShaderWithTexture, model);
}


void f_room1_window(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_bed_room1_window)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1, 0.7, 0.05));
    translate = glm::translate(model, glm::vec3(-0.8, 0.5, -4.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    shaderActivate(lightingShaderWithTexture);
    drawing_bed_room1_window.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void f_room1_bedtable(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_bed_room1_bedtable, Shader& ourShader)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.5, 0.3, 0.2));
    translate = glm::translate(model, glm::vec3(0.0, 0.0, -1.8));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    shaderActivate(lightingShaderWithTexture);
    drawing_bed_room1_bedtable.drawCubeWithTexture(lightingShaderWithTexture, model);

    shaderActivate(ourShader);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.13, 0.13, 0.13));
    translate = glm::translate(model, glm::vec3(0.17, 0.4, -1.75));
    model = alTogether * translate * scale;
    ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.5f, 0.5f, 0.5f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    
    shaderActivate(lightingShader);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.03, 0.13, 0.03));
    translate = glm::translate(model, glm::vec3(0.22, 0.3, -1.7));
    model = alTogether * translate * scale;
    drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    //shaderActivate(lightingShaderWithTexture);
    //drawing_bed_room1_bedtable.drawCubeWithTexture(lightingShaderWithTexture, model);

}

void drawMirrorWithStand(Shader& lightingShaderWithTexture, glm::mat4 alTogether, Cube& frame, Cube& mirrorSurface, Cube& stand, float tx, float ty, float tz, float rotateY, float scaleFactor) {
    // Base transformation for the mirror
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 temp = glm::mat4(1.0f);

    // Scaling matrix for the entire structure
    glm::mat4 scaleMatrix = glm::scale(glm::mat4(1.0f), glm::vec3(scaleFactor, scaleFactor, scaleFactor));
    glm::mat4 rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));

    // Adjust the Y position of the mirror and its stand
    ty += 1.0f * scaleFactor; // Bring it higher along the Y-axis, considering scaling

    glm::mat4 mirrorTransform = alTogether * glm::translate(glm::mat4(1.0f), glm::vec3(tx+0.2, ty, tz)) * rotateYMatrix * scaleMatrix;

    // Frame dimensions
    float frameWidth = 0.05f; // Thickness of the frame
    float mirrorHeight = 2.8f; // Total height of the mirror
    float mirrorWidth = 0.9f;  // Total width of the mirror
    float frameDepth = 0.05f;  // Depth of the frame

    // Mirror surface dimensions
    float surfaceHeight = mirrorHeight - frameWidth * 2;
    float surfaceWidth = mirrorWidth - frameWidth * 2;

    // Draw top frame
    temp = glm::translate(glm::mat4(1.0f), glm::vec3(-0.25f, mirrorHeight -0.9, 0.0f));
    temp = glm::scale(temp, glm::vec3(mirrorWidth, frameWidth, frameDepth));
    model = mirrorTransform * temp;
    frame.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Draw bottom frame
    temp = glm::translate(glm::mat4(1.0f), glm::vec3(-0.25f, -0.8f, 0.0f));
    temp = glm::scale(temp, glm::vec3(mirrorWidth, frameWidth, frameDepth));
    model = mirrorTransform * temp;
    frame.drawCubeWithTexture(lightingShaderWithTexture, model);


    // Draw the mirror surface
    temp = glm::translate(glm::mat4(1.0f), glm::vec3(-0.2f, -0.8f, -frameDepth / 2.0f)); // Slightly inside the frame
    temp = glm::scale(temp, glm::vec3(surfaceWidth, surfaceHeight, 0.01f)); // Thin mirror surface
    model = mirrorTransform * temp;
    mirrorSurface.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Draw the stand (angled backward)
    float standHeight = 1.6f;  // Height of the stand
    float standWidth = 0.1f;   // Thickness of the stand
    float standDepth = 0.2f;   // Depth of the stand
    float standAngle = -15.0f; // Angle of the stand

    // Stand frame
    temp = glm::translate(glm::mat4(1.0f), glm::vec3(0.1f, -1.0f, 0.4f)); // Move behind the mirror
    temp = glm::rotate(temp, glm::radians(standAngle), glm::vec3(1.0f, 0.0f, 0.0f)); // Tilt backward
    temp = glm::scale(temp, glm::vec3(standWidth, standHeight, standDepth));
    model = mirrorTransform * temp;
    stand.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Bottom connector
    float connectorHeight = 0.05f; // Height of the connector
    temp = glm::translate(glm::mat4(1.0f), glm::vec3(0.1f, -mirrorHeight / 2.0f, 0.1f)); // Bottom position
    temp = glm::scale(temp, glm::vec3(standWidth, connectorHeight, 1.0));
    model = mirrorTransform * temp;
    stand.drawCubeWithTexture(lightingShaderWithTexture, model);
}



// *************************DINING ROOM*************************************
void f_room2_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube)
{
    float baseHeight = 1.5;
    float width = 0.01;
    float length = 3.0;

    // right
    float baseHeight1 = baseHeight / 3;
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 rotateYMatrix = glm::rotate(model, glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    scale = glm::scale(model, glm::vec3(width + 0.1, baseHeight1, length + 0.1));
    translate = glm::translate(model, glm::vec3(1.6, 2*baseHeight1, -1.5));
    translate2 = glm::translate(model, glm::vec3(3.2, 0.0, -3.1));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);
    
    float baseHeight2 = baseHeight1 * 2;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width + 0.1, baseHeight2, length + 0.1 - 0.5));
    translate = glm::translate(model, glm::vec3(1.6, 0.0, -1.5));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    
    // left
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width + 0.1, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // back
    width = 3.0;
    length = 0.01;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width+0.1, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
}


void f_room2_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube)
{
    float baseHeight = 0.01;
    float width = 3.0;
    float length = 3.0;

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length + 1));
    translate = glm::translate(model, glm::vec3(0.65, 0, -1.25));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void f_room2_table(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& dining_table, Cube& mirror)
{
    float baseHeight = 0.02;
    float width = 2.2;
    float length = 0.6;

    // Table Top
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 translate4 = glm::mat4(1.0f);

    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(-1.0, 0.0, 0.1));
    translate2 = glm::translate(model, glm::vec3(3.25, 0.5, -2.25));
    model = alTogether * translate2 * translate * scale;
    shaderActivate(lightingShaderWithTexture);
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

    // Table Legs
    float legWidth = 0.05;
    float legLength = 0.5;

    // Front Left Leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(legWidth, baseHeight - 0.5, legLength));
    translate = glm::translate(model, glm::vec3(-0.5 - 0.3, 0.0, 0.2 + 0.05));
    model = alTogether * translate2 * translate * scale;
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Front Right Leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(legWidth, baseHeight - 0.5, legLength));
    translate = glm::translate(model, glm::vec3(1.1 - legWidth - 0.1 + 0.2, 0.0, 0.2 + 0.05));
    model = alTogether * translate2 * translate * scale;
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Mirror Above Table Top
    float mirrorWidth = 2.2;
    float mirrorHeight = 1.0;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(mirrorWidth, mirrorHeight, 0.1));
    translate = glm::translate(model, glm::vec3(-1.0, 0.1, 0.7));
    model = alTogether * translate2 * translate * scale;
    mirror.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Drawer Dimensions
    float drawerWidth = 0.6;
    float drawerHeight = 0.2;
    float drawerDepth = 0.5;

    // Drawer 1 (Left)
    // Drawer 1 Bottom
    translate3 = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 0.0, bathroom_door_translate));
    translate4 = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 0.0, bathroom_door_translate*25.0));

    model = glm::mat4(1.0f);
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(drawerWidth, 0.02, drawerDepth));
    translate = glm::translate(glm::mat4(1.0f), glm::vec3(-0.51, -0.22, 0.2));
    model = alTogether * translate2 * translate * scale;
    model = model * translate3;
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Drawer 1 Left Wall
    model = glm::mat4(1.0f);
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.02, drawerHeight, drawerDepth));
    translate = glm::translate(glm::mat4(1.0f), glm::vec3(-0.2 - drawerWidth / 2 + 0.01, -0.22, 0.2));
    model = alTogether * translate2 * translate * scale;
    model = model * translate3;
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Drawer 1 Right Wall
    model = glm::mat4(1.0f);
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.02, drawerHeight, drawerDepth));
    translate = glm::translate(glm::mat4(1.0f), glm::vec3(-0.2 + drawerWidth / 2 - 0.01, -0.22, 0.2));
    model = alTogether * translate2 * translate * scale;
    model = model * translate3;
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Drawer 1 Back Wall
    model = glm::mat4(1.0f);
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(drawerWidth, drawerHeight, 0.02));
    translate = glm::translate(glm::mat4(1.0f), glm::vec3(-0.5, -0.22, 0.2));
    model = alTogether * translate2 * translate * scale;
    model = model * translate4;
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);
    // Drawer 1 Back Wall
    model = glm::mat4(1.0f);
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(drawerWidth, drawerHeight, 0.02));
    translate = glm::translate(glm::mat4(1.0f), glm::vec3(-0.5, -0.22, 0.7));
    model = alTogether * translate2 * translate * scale;
    model = model * translate4;
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);
  
    // Drawer 2 (Right)
    // Drawer 2 Bottom
    model = glm::mat4(1.0f);
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(drawerWidth, 0.02, drawerDepth));
    translate = glm::translate(glm::mat4(1.0f), glm::vec3(0.4, -0.22, 0.2));
    model = alTogether * translate2 * translate * scale;
    model = model * translate3;
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);
    
    // Drawer 2 Left Wall
    model = glm::mat4(1.0f);
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.02, drawerHeight, drawerDepth));
    translate = glm::translate(glm::mat4(1.0f), glm::vec3(0.7 - drawerWidth / 2 + 0.01, -0.22, 0.2));
    model = alTogether * translate2 * translate * scale;
    model = model * translate3;
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Drawer 2 Right Wall
    model = glm::mat4(1.0f);
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.02, drawerHeight, drawerDepth));
    translate = glm::translate(glm::mat4(1.0f), glm::vec3(0.7 + drawerWidth / 2 - 0.01, -0.22, 0.2));
    model = alTogether * translate2 * translate * scale;
    model = model * translate3;
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Drawer 2 Back Wall
    model = glm::mat4(1.0f);
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(drawerWidth, drawerHeight, 0.02));
    translate = glm::translate(glm::mat4(1.0f), glm::vec3(0.4, -0.22, 0.2));
    model = alTogether * translate2 * translate * scale;
    model = model * translate4;
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Drawer 2 Back Wall
    model = glm::mat4(1.0f);
    scale = glm::scale(glm::mat4(1.0f), glm::vec3(drawerWidth, drawerHeight, 0.02));
    translate = glm::translate(glm::mat4(1.0f), glm::vec3(0.4, -0.22, 0.7));
    model = alTogether * translate2 * translate * scale;
    model = model * translate4;
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);


}

void f_room2_table2(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& dining_table, Cube& mirror)
{
    float baseHeight = 0.02;
    float width = 3.0;
    float length = 0.6;

    // Table Top
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(length, baseHeight, width));
    //translate = glm::translate(model, glm::vec3(1.0, 0.0, -0.1));
    translate2 = glm::translate(model, glm::vec3(-1.2, 0.65, -4.5));
    model = alTogether * translate2 * scale;
    shaderActivate(lightingShaderWithTexture);
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

    // Table Legs
    float legWidth = 0.05;
    float legLength = 0.65;

    // Front Left Leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(legLength, baseHeight , legWidth));
    translate = glm::translate(model, glm::vec3(0.2, 0.0, 0.45));
    model = alTogether * translate2 * translate * scale;
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Front Right Leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(legLength, baseHeight , legWidth));
    translate = glm::translate(model, glm::vec3(0.2, 0.0, 2.8));
    model = alTogether * translate2 * translate * scale;
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Mirror Above Table Top
    float mirrorWidth = 2.8;
    float mirrorHeight = 1.0;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, mirrorHeight, mirrorWidth));
    translate = glm::translate(model, glm::vec3(0.0, 0.1, 0.2));
    model = alTogether * translate2 * translate * scale;
    mirror.drawCubeWithTexture(lightingShaderWithTexture, model);

}

void f_room2_frize(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_dining_frize)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.5, 1.0, 0.3));
    translate = glm::translate(model, glm::vec3(1.75, 0.0, -4.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    shaderActivate(lightingShaderWithTexture);
    drawing_dining_frize.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void f_room2_window(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& kitchen_window)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.02, 0.6, 0.9));
    translate = glm::translate(model, glm::vec3(4.55, 0.6, -2.75));
    model = alTogether * translate * scale;
    //kitchen_window.drawCubeWithTexture(lightingShaderWithTexture, model);
}


glm::mat4 transform(float tx, float ty, float tz, float sx, float sy, float sz) {
    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;
    translateMatrix = glm::translate(identityMatrix, glm::vec3(tx, ty, tz));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(sx, sy, sz));
    model = translateMatrix * scaleMatrix;
    return model;
}


void chair(unsigned int& cubeVAO, Shader& lightingShaderWithTexture, glm::mat4 alTogether, Cube& chair_texture, Cube& dining_frize, float tx, float ty, float tz, float rotateY)
{
    // Chair components' dimensions (relative to its local space)
    glm::vec3 seatPos(0.0f, 0.0f, 0.0f), seatDim(1.3f, 0.4f, 1.3f);
    glm::vec3 backrestPos(0.0f, 0.1f, -0.2f), backrestDim(1.3f, 1.0f, 0.2f);
    glm::vec3 leftArmPos(-0.2f, 0.1f, 0.0f), leftArmDim(0.2f, 0.4f, 1.0f);
    glm::vec3 rightArmPos(1.1f, 0.1f, 0.0f), rightArmDim(0.2f, 0.4f, 1.0f);
    glm::vec3 polePos(0.6f, -1.1f, 0.0f), poleDim(0.2f, 1.3f, 0.2f);
    glm::vec3 basePos(0.0f, -1.2f, 0.0f), baseDim(1.5f, 0.1f, 1.5f);

    // Calculate the center of the seat
    glm::vec3 centerOfSeat = seatPos + glm::vec3(seatDim.x / 2.0f, seatDim.y / 2.0f, seatDim.z / 2.0f);

    // Transformation matrices
    glm::mat4 scaleMatrix = glm::scale(glm::mat4(1.0f), glm::vec3(0.3f, 0.3f, 0.3f));
    glm::mat4 toCenterOfSeat = glm::translate(glm::mat4(1.0f), -centerOfSeat);
    glm::mat4 rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(rotateY), glm::vec3(0.0f, 1.0f, 0.0f));
    glm::mat4 fromCenterOfSeat = glm::translate(glm::mat4(1.0f), centerOfSeat);
    glm::mat4 centerSeatRotation = fromCenterOfSeat * rotateYMatrix * toCenterOfSeat;

    // Combine with global transformations
    glm::mat4 chairTransform = alTogether * glm::translate(glm::mat4(1.0f), glm::vec3(tx, ty, tz)) * scaleMatrix * centerSeatRotation;

    // Seat
    glm::mat4 temp = transform(seatPos.x, seatPos.y, seatPos.z, seatDim.x, seatDim.y, seatDim.z);
    glm::mat4 model = chairTransform * temp;
    chair_texture.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Backrest
    temp = transform(backrestPos.x, backrestPos.y, backrestPos.z, backrestDim.x, backrestDim.y, backrestDim.z);
    model = chairTransform * temp;
    chair_texture.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Left Armrest
    temp = transform(leftArmPos.x, leftArmPos.y, leftArmPos.z, leftArmDim.x, leftArmDim.y, leftArmDim.z);
    model = chairTransform * temp;
    chair_texture.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Right Armrest
    temp = transform(rightArmPos.x, rightArmPos.y, rightArmPos.z, rightArmDim.x, rightArmDim.y, rightArmDim.z);
    model = chairTransform * temp;
    chair_texture.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Hydraulic Support (Pole)
    temp = transform(polePos.x, polePos.y, polePos.z, poleDim.x, poleDim.y, poleDim.z);
    model = chairTransform * temp;
    dining_frize.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Base
    temp = transform(basePos.x, basePos.y, basePos.z, baseDim.x, baseDim.y, baseDim.z);
    model = chairTransform * temp;
    dining_frize.drawCubeWithTexture(lightingShaderWithTexture, model);
}

void chair2(unsigned int& cubeVAO, Shader& lightingShaderWithTexture, glm::mat4 alTogether, Cube& chair_texture, Cube& dining_frize, float tx, float ty, float tz, float rotateY)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 temp = glm::mat4(1.0f);
    glm::mat4 scaleMatrix = glm::scale(glm::mat4(1.0f), glm::vec3(0.3, 0.3, 0.3));
    glm::mat4 rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(rotateY), glm::vec3(0.0f, 1.0f, 0.0f));

    // Seat
    temp = transform(tx+0.9, ty, tz+0.5, 1.3, 0.25, 1.7); // Wider and slightly flatter
    model = alTogether * rotateYMatrix * scaleMatrix * temp;
   // chair_texture.drawCubeWithTexture(lightingShaderWithTexture, model);
    seat->draw(lightingShaderWithTexture, model);

    // Backrest (Dynamic Tilt)
    temp = transform(tx + 0.9, ty + 0.8, tz - 0.5, 1.2, 1.3, 0.15);
    temp = glm::rotate(temp, glm::radians(backrestTiltAngle), glm::vec3(1.0f, 0.0f, 0.0f)); // Tilt backward dynamically
    model = alTogether * rotateYMatrix * scaleMatrix * temp;
    backrest2->draw(lightingShaderWithTexture, model);

    // Left Armrest
    temp = transform(tx+0.1 , ty + 0.4, tz+0.2, 0.2, 0.4, 1.4); // Adjust position and size
    model = alTogether * rotateYMatrix * scaleMatrix * temp;
    //chair_texture.drawCubeWithTexture(lightingShaderWithTexture, model);
    backrest2->draw(lightingShaderWithTexture, model);

    // Right Armrest
    temp = transform(tx + 1.7, ty + 0.4, tz+0.2, 0.2, 0.4, 1.4); // Adjust position and size
    model = alTogether * rotateYMatrix * scaleMatrix * temp;
    //chair_texture.drawCubeWithTexture(lightingShaderWithTexture, model);
    backrest2->draw(lightingShaderWithTexture, model);

    // Hydraulic Support (Pole)
    temp = transform(tx +0.8, ty-1.3 , tz+0.4, 0.2, 1.2, 0.2); // Centered pole
    model = alTogether * rotateYMatrix * scaleMatrix * temp;
    dining_frize.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Base
    temp = transform(tx, ty - 1.5, tz, 2.0, 0.2, 2.0); // Larger and flatter base
    model = alTogether * rotateYMatrix * scaleMatrix * temp;
    dining_frize.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Footrest
    temp = transform(tx+0.4, ty - 0.8, tz + 1.7, 1.2, 0.1, 0.5); // Footrest below seat
    model = alTogether * rotateYMatrix * scaleMatrix * temp;
    dining_frize.drawCubeWithTexture(lightingShaderWithTexture, model);

    // Handle for the Footrest
    temp = transform(tx + 0.9, ty - 0.8, tz + 0.4, 0.1 ,0.1, 1.3); // Handle connecting seat and footrest
    model = alTogether * rotateYMatrix * scaleMatrix * temp;
    dining_frize.drawCubeWithTexture(lightingShaderWithTexture, model);

}

void f_room2_light(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Shader& ourShader)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    // point light
    shaderActivate(ourShader);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.15, 0.15));
    translate = glm::translate(model, glm::vec3(1.6f, 1.3f, -3.1f));
    model = alTogether * translate * scale;
    ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    shaderActivate(lightingShader);
}
// **********************************       Kitchen           ***********************************
void f_kitchen_surface(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& kitchen_surface, Cube& kitchen_surface_top, Cube& kitchen_cupboard1, Cube& kitchen_cupboard2, Cube& kitchen_back_texture)
{
    // ground
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.0, 0.5, 0.4));
    translate = glm::translate(model, glm::vec3(2.75, 0.0, -4.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    shaderActivate(lightingShaderWithTexture);
    kitchen_surface.drawCubeWithTexture(lightingShaderWithTexture, model);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.0, 0.02, 0.4));
    translate = glm::translate(model, glm::vec3(2.75, 0.5, -4.5));
    model = alTogether * translate * scale;
    kitchen_surface_top.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.4, 0.5, 0.9));
    translate = glm::translate(model, glm::vec3(4.25, 0.0, -4.1));
    model = alTogether * translate * scale;
    kitchen_surface.drawCubeWithTexture(lightingShaderWithTexture, model);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.4, 0.02, 0.9));
    translate = glm::translate(model, glm::vec3(4.25, 0.5, -4.1));
    model = alTogether * translate * scale;
    kitchen_surface_top.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.0, 0.5, 0.4));
    translate = glm::translate(model, glm::vec3(2.25, 0.0, -3.6));
    model = alTogether * translate * scale;
    kitchen_surface.drawCubeWithTexture(lightingShaderWithTexture, model);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.0, 0.02, 0.4));
    translate = glm::translate(model, glm::vec3(2.25, 0.5, -3.6));
    model = alTogether * translate * scale;
    kitchen_surface_top.drawCubeWithTexture(lightingShaderWithTexture, model);

    // cupboard
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.0, 0.5, 0.4));
    translate = glm::translate(model, glm::vec3(2.75, 1.0, -4.5));
    model = alTogether * translate * scale;
    kitchen_cupboard1.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.4, 0.5, 0.9));
    translate = glm::translate(model, glm::vec3(4.25, 1.0, -4.1));
    model = alTogether * translate * scale;
    kitchen_cupboard2.drawCubeWithTexture(lightingShaderWithTexture, model);

    // kitchen_back
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.5, 0.4, 0.02));
    translate = glm::translate(model, glm::vec3(2.75, 0.5, -4.5));
    model = alTogether * translate * scale;
    kitchen_back_texture.drawCubeWithTexture(lightingShaderWithTexture, model);


    shaderActivate(lightingShader);
}


void f_kitchen_stove(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& kitchen_surface)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.6, 0.06, 0.4));
    translate = glm::translate(model, glm::vec3(3.3, 0.52, -3.6));
    model = alTogether * translate * scale;
    kitchen_surface.drawCubeWithTexture(lightingShaderWithTexture, model);
    
}
// **********************************       BED ROOM 2        ************************************
void f_room3_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube)
{
    float baseHeight = 1.5;
    float width = 0.01;
    float length = 3.0;

    // right
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width + 0.1, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(1.6, 0.0, -1.5));
    translate2 = glm::translate(model, glm::vec3(3.0, 0.0, 0.0));
    model = alTogether * translate2 * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);


    // back
    // back top
    float baseHeight1 = baseHeight / 3;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(length + 0.1, baseHeight1, width + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 2 * baseHeight1, 1.5));
    glm::mat4 translate3 = glm::translate(model, glm::vec3(0.0, 0.0, -1.75));
    model = alTogether * translate3 * translate2 * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // front bottom
    float baseHeight2 = baseHeight1 * 2;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(length - 0.5, baseHeight2, width + 0.1));
    translate = glm::translate(model, glm::vec3(-0.9, 0.0, 1.5));
    model = alTogether * translate3 * translate2 * translate * scale;/*
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);*/

    // back
    width = 3.0;
    length = 0.01;
    //model = glm::mat4(1.0f);
    //scale = glm::scale(model, glm::vec3(width, baseHeight, length + 0.1));
    //translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    //model = alTogether * translate * scale;
    ////drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //

    // front
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width + 0.1, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 0.0, 1.5));
    model = alTogether * translate2 * translate * scale;
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    

    shaderActivate(lightingShader);

}


void f_room3_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube)
{
    float baseHeight = 0.01;
    float width = 3.0;
    float length = 3.0;

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length + 1));
    translate = glm::translate(model, glm::vec3(0.65, 0, -0.5));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void f_room3_book(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_dining_frize)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.8, 0.6, 0.2));
    translate = glm::translate(model, glm::vec3(1.75, 0.5, 1.3));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    shaderActivate(lightingShaderWithTexture);
    drawing_dining_frize.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void f_room3_light(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Shader& ourShader)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    // point light
    shaderActivate(ourShader);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.15, 0.15));
    translate = glm::translate(model, glm::vec3(1.6f, 1.3f, 0.5f));
    model = alTogether * translate * scale;
    ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    shaderActivate(lightingShader);
}


void f_room3_bed(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& bed_sheet, Cube& bed_pillow, Cube& bed_texture, Cube& blanket_texture)
{
    float baseHeight = 0.3;
    float width = 0.7;
    float length = 2;
    float pillowWidth = 0.25;
    float pillowLength = 0.2;
    float blanketWidth = 0.7;
    float blanketLength = 0.5;
    float headHeight = 0.6;

    //base
    glm::mat4 translate_origin = glm::mat4(1.0f);
    glm::mat4 translate_location = glm::mat4(1.0f);
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 rotateM = glm::rotate(model, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    translate_origin = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    translate_location = glm::translate(model, glm::vec3(6.3, 0.0, -0.45));
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    translate3 = glm::translate(model, glm::vec3(-0.75, 0, -2.7));
    model = alTogether * translate_location * rotateM * translate3 * scale * translate;
    shaderActivate(lightingShaderWithTexture);
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    bed_texture.drawCubeWithTexture(lightingShaderWithTexture, model);

    //foam
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, baseHeight, 0));
    scale = glm::scale(model, glm::vec3(width, 0.06, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate_location * rotateM * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.804, 0.361, 0.361);
    bed_sheet.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

    //pillow 1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3((width / 2) - (0.1 + pillowWidth / 2), baseHeight + 1 * 0.06, (length / 2) - (0.025 + pillowWidth / 2)));
    scale = glm::scale(model, glm::vec3(pillowWidth, 0.04, pillowLength));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate_location * rotateM * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 1, 0.647, 0);
    shaderActivate(lightingShaderWithTexture);
    bed_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);

    //pillow 2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3((-width / 2) + (0.1 + pillowWidth / 2), baseHeight + 1 * 0.06, (length / 2) - (0.025 + pillowWidth / 2)));
    scale = glm::scale(model, glm::vec3(pillowWidth, 0.04, pillowLength));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate_location * rotateM * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 1, 0.647, 0);
    bed_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

    //blanket
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, baseHeight + 1 * 0.06, -(length / 2 - 0.025) + blanketLength / 2));
    scale = glm::scale(model, glm::vec3(blanketWidth, 0.015, blanketLength));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate_location * rotateM * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.541, 0.169, 0.886);
    blanket_texture.drawCubeWithTexture(lightingShaderWithTexture, model);

    //head
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, 0, (length / 2 - 0.02 / 2) + 0.02));
    scale = glm::scale(model, glm::vec3(width, headHeight, 0.06));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate_location * rotateM * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    bed_texture.drawCubeWithTexture(lightingShaderWithTexture, model);
}


void f_room3_dressing_table(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& dressing_right, Cube& dressing_bottom, Cube& dressing_mirror)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    // mirror
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.6, 0.3));
    translate = glm::translate(model, glm::vec3(4.55f, 0.35f, 0.9f));
    model = alTogether * translate * scale;
    dressing_mirror.drawCubeWithTexture(lightingShaderWithTexture, model);


}
// *********************************   ROOF    *****************************************************
void f_roof(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_roof)
{
    float baseHeight = 0.01;
    float width = 3.0;
    float length = 3.0;

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2*(width + 1), baseHeight, 2*(length + 1)));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    translate2 = glm::translate(model, glm::vec3(1.5, 1.5, -1.5));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    drawing_roof.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}

// *************************BATHROOM****************************
void f_bathroom_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& bathroom_top, Cube& bathroom_door, Cube& bathroom_tiles)
{

    float baseHeight = 1.5;
    float width = 0.01;
    float length = 3.0;
    float base_top = baseHeight / 3.0;
    shaderActivate(lightingShaderWithTexture);

    // top
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, base_top, length + 0.1 - 1.75));
    translate = glm::translate(model, glm::vec3(2.5, 2*base_top, -1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    bathroom_top.drawCubeWithTexture(lightingShaderWithTexture, model);

    // left_door
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, 2*base_top, (length + 0.1 - 1.75)/2.0));
    translate = glm::translate(model, glm::vec3(2.5, 0.0, -1.5 + bathroom_door_translate));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    bathroom_door.drawCubeWithTexture(lightingShaderWithTexture, model);

    
    // right_door
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, 2 * base_top, (length + 0.1 - 1.75) / 2.0));
    translate = glm::translate(model, glm::vec3(2.5, 0.0, -1.5+ (length + 0.1 - 1.75) / 2.0));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    bathroom_door.drawCubeWithTexture(lightingShaderWithTexture, model);

    // back_wall
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, 2 * base_top, length + 0.1 - 1.75));
    translate = glm::translate(model, glm::vec3(2.5 + 2.0, 0.0, -1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    bathroom_tiles.drawCubeWithTexture(lightingShaderWithTexture, model);

    // left_wall
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(length + 0.1 - 1.25, 2 * base_top, width));
    translate = glm::translate(model, glm::vec3(2.6, 0.0, -1.5+0.15));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    bathroom_tiles.drawCubeWithTexture(lightingShaderWithTexture, model);

    // right_wall
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(length + 0.1 - 1.25, 2 * base_top, width));
    translate = glm::translate(model, glm::vec3(2.6, 0.0, -0.3));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    bathroom_tiles.drawCubeWithTexture(lightingShaderWithTexture, model);

    // top_wall
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.0, width, 2*max_bathroom_door_translate));
    translate = glm::translate(model, glm::vec3(2.6, 2*base_top, -1.65 + 0.15));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    bathroom_tiles.drawCubeWithTexture(lightingShaderWithTexture, model);

    // bottom wall
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.0, width, 2 * max_bathroom_door_translate));
    translate = glm::translate(model, glm::vec3(2.6, 0.1, -1.65 + 0.15));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    bathroom_tiles.drawCubeWithTexture(lightingShaderWithTexture, model);


    shaderActivate(lightingShader);
}

void f_bathroom_light(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Shader& ourShader)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    // point light
    shaderActivate(ourShader);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.15, 0.15));
    translate = glm::translate(model, glm::vec3(2.5 + 1.9, 0.8, -0.9));
    model = alTogether * translate * scale;
    ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    shaderActivate(lightingShader);
}

void f_toilet(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& toilet)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    // toilet
    shaderActivate(lightingShaderWithTexture);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.4, 0.05, 0.6));
    translate = glm::translate(model, glm::vec3(2.5 + 1.5, 0.1, -0.95));
    model = alTogether * translate * scale;
    /*ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);*/
    toilet.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}

void f_bathroom_bucket(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& tex)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 translate4 = glm::mat4(1.0f);
    glm::mat4 translate5 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 scale2 = glm::mat4(1.0f);
    glm::mat4 scale3 = glm::mat4(1.0f);

    // toilet
    shaderActivate(lightingShaderWithTexture);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.20, 0.25, 0.35));
    scale2 = glm::scale(model, glm::vec3(0.20, 0.15, 0.25));
    scale3 = glm::scale(model, glm::vec3(0.05, 0.25, 0.05));
    translate = glm::translate(model, glm::vec3(2.5 + 1.0, 0.0, -1.3));
    translate2 = glm::translate(model, glm::vec3(-0.2, -0.8, 0.0));
    translate3 = glm::translate(model, glm::vec3(0.0, 0.08, 0.6));
    translate4 = glm::translate(model, glm::vec3(0.1, 0.0, 0.6));
    translate5 = glm::translate(model, glm::vec3(-0.15, 0.0, 0.6));
    glm::mat4 rotateM = glm::rotate(model, glm::radians(180.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    glm::mat4 rotateN = glm::rotate(model, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    glm::mat4 rotateo = glm::rotate(model, glm::radians(70.0f), glm::vec3(0.0f, 0.0f, 1.0f));

    model = alTogether * translate * rotateM * rotateN* rotateo* translate2 * scale;
    
    /*ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);*/

  
    

    tex.drawCubeWithTexture(lightingShaderWithTexture, translate * translate4   *scale3);
    tex.drawCubeWithTexture(lightingShaderWithTexture, translate  *translate5 * scale3);
    bucket->draw(lightingShaderWithTexture, model);
    backrest->draw(lightingShaderWithTexture, translate * translate3* rotateN * scale2);
    //shaderActivate(lightingShader);

   
   
}


void getCurrentTime(int& hours, int& minutes, int& seconds) {
    time_t currentTime = time(nullptr); // Get current UNIX timestamp
    struct tm* timeinfo;
    timeinfo = localtime(&currentTime);

    seconds = timeinfo->tm_sec;
    minutes = timeinfo->tm_min;
    hours = timeinfo->tm_hour;
}
// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        camera.ProcessKeyboard(FORWARD, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        camera.ProcessKeyboard(LEFT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        camera.ProcessKeyboard(RIGHT, deltaTime);
    }

    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) {
        camera.ProcessKeyboard(UP, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
        camera.ProcessKeyboard(DOWN, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_T) == GLFW_PRESS) {
        camera.ProcessKeyboard(P_UP, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS) {
        camera.ProcessKeyboard(P_DOWN, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS) {
        camera.ProcessKeyboard(Y_LEFT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_H) == GLFW_PRESS) {
        camera.ProcessKeyboard(Y_RIGHT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) {
        camera.ProcessKeyboard(R_LEFT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_J) == GLFW_PRESS) {
        camera.ProcessKeyboard(R_RIGHT, deltaTime);
    }

    /*if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
    {
        if (rotateAxis_X) rotateAngle_X -= 0.1;
        else if (rotateAxis_Y) rotateAngle_Y -= 0.1;
        else rotateAngle_Z -= 0.1;
    }*/
    /*if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) translate_Y += 0.001;*/
    //if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS) translate_Y -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS) translate_X += 0.001;
    /*if (glfwGetKey(window, GLFW_KEY_J) == GLFW_PRESS) translate_X -= 0.001;*/
    if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) translate_Z += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) translate_Z -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS) scale_X += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS) scale_X -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS) scale_Y += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS) scale_Y -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS) scale_Z += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS) scale_Z -= 0.001;

    if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS)
    {
        rotateAngle_X += 0.1;
        rotateAxis_X = 1.0;
        rotateAxis_Y = 0.0;
        rotateAxis_Z = 0.0;
    }
    if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS)
    {
        rotateAngle_Y += 0.1;
        rotateAxis_X = 0.0;
        rotateAxis_Y = 1.0;
        rotateAxis_Z = 0.0;
    }
    if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS)
    {
        rotateAngle_Z += 0.1;
        rotateAxis_X = 0.0;
        rotateAxis_Y = 0.0;
        rotateAxis_Z = 1.0;
    }

    /*if (glfwGetKey(window, GLFW_KEY_H) == GLFW_PRESS)
    {
        eyeX += 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }*/
    /*if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS)
    {
        eyeX -= 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }*/
    /*if (glfwGetKey(window, GLFW_KEY_T) == GLFW_PRESS)
    {
        eyeZ += 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }*/
    if (glfwGetKey(window, GLFW_KEY_G) == GLFW_PRESS)
    {
        eyeZ -= 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
    {
        eyeY += 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }
    /*if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
    {
        eyeY -= 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }*/
    // front door
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
    {
        if (bathroom_door_translate == 0.0f)
        {
            f_door += 1;
            f_door = min(70.0f, f_door);
        }
        else
        {
            spongeup += 0.4;
        }
        
    }
    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
    {
        
        if (bathroom_door_translate == 0.0f)
        {
            f_door -= 1;
            f_door = max(0.0f, f_door);
        }
        else
        {
            spongeup -= 0.4;
        }
    }
    // back door
    if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS)
    {
        b1_door += 1;
        b1_door = min(80.0f, b1_door);
    }
    if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS)
    {
        b1_door -= 1;
        b1_door = max(0.0f, b1_door);
    }
    // bathroom door
    if (glfwGetKey(window, GLFW_KEY_5) == GLFW_PRESS)
    {
        bathroom_door_translate -= 0.02;
        bathroom_door_translate = min(max_bathroom_door_translate, bathroom_door_translate);
    }
    if (glfwGetKey(window, GLFW_KEY_6) == GLFW_PRESS)
    {
        bathroom_door_translate += 0.01;
        bathroom_door_translate = max(0.0f, bathroom_door_translate);
    }
    //chair
    if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS) {
        chairRotation += 1.0f; // Toggle rotation on/off
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT_BRACKET) == GLFW_PRESS) { // Backtick key `\``
        
       if(backrestTiltAngle>=-80.0f) backrestTiltAngle -= 0.5f; // Toggle tilting
    }
    if (glfwGetKey(window, GLFW_KEY_RIGHT_BRACKET) == GLFW_PRESS) { // Backtick key `\``
        backrestTiltAngle = -10.0f; // Toggle tilting
    }
    
    
    // fan
    if (glfwGetKey(window, GLFW_KEY_0) == GLFW_PRESS)
    {
        isFanOn = !isFanOn; // Toggle the switch state
        if(isFanOn) rotateFan += 5.0f;
        
    }
    
    if (glfwGetKey(window, GLFW_KEY_7) == GLFW_PRESS)
    {
        pointlight1.turnAmbientOn();
        pointlight3.turnAmbientOn();
        drawing_light.turnAmbientOn();
        bed_room1_light.turnAmbientOn();
        dining_light.turnAmbientOn();
        bed_room2_light.turnAmbientOn();
        bathroom_light.turnAmbientOn();

        pointlight1.turnDiffuseOff();
        pointlight3.turnDiffuseOff();
        drawing_light.turnDiffuseOff();
        bed_room1_light.turnDiffuseOff();
        dining_light.turnDiffuseOff();
        bed_room2_light.turnDiffuseOff();
        bathroom_light.turnDiffuseOff();

        pointlight1.turnSpecularOff();
        pointlight3.turnSpecularOff();
        drawing_light.turnSpecularOff();
        bed_room1_light.turnSpecularOff();
        dining_light.turnSpecularOff();
        bed_room2_light.turnSpecularOff();
        bathroom_light.turnSpecularOff();
    }
    
    if (glfwGetKey(window, GLFW_KEY_8) == GLFW_PRESS)
    {
        pointlight1.turnDiffuseOn();
        pointlight3.turnDiffuseOn();
        drawing_light.turnDiffuseOn();
        bed_room1_light.turnDiffuseOn();
        dining_light.turnDiffuseOn();
        bed_room2_light.turnDiffuseOn();
        bathroom_light.turnDiffuseOn();

        pointlight1.turnAmbientOff();
        pointlight3.turnAmbientOff();
        drawing_light.turnAmbientOff();
        bed_room1_light.turnAmbientOff();
        dining_light.turnAmbientOff();
        bed_room2_light.turnAmbientOff();
        bathroom_light.turnAmbientOff();

        pointlight1.turnSpecularOff();
        pointlight3.turnSpecularOff();
        drawing_light.turnSpecularOff();
        bed_room1_light.turnSpecularOff();
        dining_light.turnSpecularOff();
        bed_room2_light.turnSpecularOff();
        bathroom_light.turnSpecularOff();
    }

    if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS)
    {
        pointlight1.turnDiffuseOff();
        pointlight3.turnDiffuseOff();
        drawing_light.turnDiffuseOff();
        bed_room1_light.turnDiffuseOff();
        dining_light.turnDiffuseOff();
        bed_room2_light.turnDiffuseOff();
        bathroom_light.turnDiffuseOff();

        pointlight1.turnAmbientOff();
        pointlight3.turnAmbientOff();
        drawing_light.turnAmbientOff();
        bed_room1_light.turnAmbientOff();
        dining_light.turnAmbientOff();
        bed_room2_light.turnAmbientOff();
        bathroom_light.turnAmbientOff();

        pointlight1.turnSpecularOn();
        pointlight3.turnSpecularOn();
        drawing_light.turnSpecularOn();
        bed_room1_light.turnSpecularOn();
        dining_light.turnSpecularOn();
        bed_room2_light.turnSpecularOn();
        bathroom_light.turnSpecularOn();
    }

    if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS)
    {
        pointlight1.turnDiffuseOn();
        pointlight3.turnDiffuseOn();
        drawing_light.turnDiffuseOn();
        bed_room1_light.turnDiffuseOn();
        dining_light.turnDiffuseOn();
        bed_room2_light.turnDiffuseOn();
        bathroom_light.turnDiffuseOn();

        pointlight1.turnAmbientOn();
        pointlight3.turnAmbientOn();
        drawing_light.turnAmbientOn();
        bed_room1_light.turnAmbientOn();
        dining_light.turnAmbientOn();
        bed_room2_light.turnAmbientOn();
        bathroom_light.turnAmbientOn();

        pointlight1.turnSpecularOn();
        pointlight3.turnSpecularOn();
        drawing_light.turnSpecularOn();
        bed_room1_light.turnSpecularOn();
        dining_light.turnSpecularOn();
        bed_room2_light.turnSpecularOn();
        bathroom_light.turnSpecularOn();
    }


   // if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
   // {
   //     if (onOffPointToggle)
   //     {
   //         pointlight1.turnOff();
   //         
   //         onOffPointToggle = false;
   //     }
   //     else
   //     {
   //         pointlight1.turnOn();
   //       
   //         onOffPointToggle = true;
   //     }
   //    // pointlight3.turnOff();
   //    // pointlight4.turnOff();

   // }
   // 

   // if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
   // {
   //     
   //     if (onOffSpotToggle)
   //     {
   //        
   //         pointlight2.turnOff();
   //         onOffSpotToggle = false;
   //     }
   //     else
   //     {
   //         pointlight2.turnOn();
   //         onOffSpotToggle = true;
   //     }
   // }

   // if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS)
   // {

   //     if (onOffDirectToggle)
   //     {

   //         pointlight3.turnOff();
   //         onOffDirectToggle = false;
   //     }
   //     else
   //     {
   //         pointlight3.turnOn();
   //         onOffDirectToggle = true;
   //     }
   // }
   // 
   // if (glfwGetKey(window, GLFW_KEY_5) == GLFW_PRESS)
   // {
   //     pointlight1.turnAmbientOn();
   //     pointlight2.turnAmbientOn();
   //    // pointlight3.turnAmbientOn();
   //    // pointlight4.turnAmbientOn();
   // }
   // if (glfwGetKey(window, GLFW_KEY_6) == GLFW_PRESS)
   // {
   //     pointlight1.turnAmbientOff();
   //     pointlight2.turnAmbientOff();
   //   //  pointlight3.turnAmbientOff();
   //   //  pointlight4.turnAmbientOff();
   // }
   // if (glfwGetKey(window, GLFW_KEY_7) == GLFW_PRESS)
   // {
   //     pointlight1.turnDiffuseOn();
   //     pointlight2.turnDiffuseOn();
   //  //   pointlight3.turnDiffuseOn();
   // //    pointlight4.turnDiffuseOn();
   // }
   // if (glfwGetKey(window, GLFW_KEY_8) == GLFW_PRESS)
   // {
   //     pointlight1.turnDiffuseOff();
   //     pointlight2.turnDiffuseOff();
   ////     pointlight3.turnDiffuseOff();
   // //    pointlight4.turnDiffuseOff();
   // }
   // if (glfwGetKey(window, GLFW_KEY_9) == GLFW_PRESS)
   // {
   //     pointlight1.turnSpecularOn();
   //     pointlight2.turnSpecularOn();
   // //    pointlight3.turnSpecularOn();
   // //    pointlight4.turnSpecularOn();
   // }
   // if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS)
   // {
   //     pointlight1.turnSpecularOff();
   //     pointlight2.turnSpecularOff();
   ////     pointlight3.turnSpecularOff();
   // //    pointlight4.turnDiffuseOff();
   // }
    if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS)
    {
        //pointlight2.turnOn();
        drawing_light.turnOn();
        bed_room1_light.turnOn();
        dining_light.turnOn();
        bed_room2_light.turnOn();
        bathroom_light.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS)
    {
        //pointlight2.turnOff();
        drawing_light.turnOff();
        bed_room1_light.turnOff();
        dining_light.turnOff();
        bed_room2_light.turnOff();
        bathroom_light.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS)
    {
        pointlight3.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS)
    {
        pointlight3.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS)
    {
        pointlight1.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS)
    {
        pointlight1.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    //if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
    //{
    //    /*if (diffuseToggle)
    //    {*/
    //    if (pointlight1.isOn())
    //        pointlight1.turnAmbientOn();
    //    if (pointlight2.isOn())
    //        pointlight2.turnAmbientOn();
    //    if (pointlight3.isOn())
    //        pointlight3.turnAmbientOn();
    //    //pointlight4.turnDiffuseOn();
    //    //diffuseToggle = !diffuseToggle;
    ////}
    //}
    //if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
    //{
    //    /*if (diffuseToggle)
    //    {*/
    //    if (pointlight1.isOn())
    //        pointlight1.turnAmbientOff();
    //    if (pointlight2.isOn())
    //        pointlight2.turnAmbientOff();
    //    if (pointlight3.isOn())
    //        pointlight3.turnAmbientOff();
    //    //pointlight4.turnDiffuseOff();
    //    //diffuseToggle = !diffuseToggle;
    ////}
    //}

    //if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS)
    //{
    //    /*if (diffuseToggle)
    //    {*/
    //    if (pointlight1.isOn())
    //        pointlight1.turnDiffuseOn();
    //    if (pointlight2.isOn())
    //        pointlight2.turnDiffuseOn();
    //    if (pointlight3.isOn())
    //        pointlight3.turnDiffuseOn();
    //    //pointlight4.turnAmbientOn();
    //    //diffuseToggle = !diffuseToggle;
    //    //}
    //}
    //if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS)
    //{
    //    /*if (diffuseToggle)
    //    {*/
    //    if (pointlight1.isOn())
    //        pointlight1.turnDiffuseOff();
    //    if (pointlight2.isOn())
    //        pointlight2.turnDiffuseOff();
    //    if (pointlight3.isOn())
    //        pointlight3.turnDiffuseOff();
    //    //diffuseToggle = !diffuseToggle;
    //    //}
    //}


    if (glfwGetKey(window, GLFW_KEY_5) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        if (pointlight1.isOn())
            pointlight1.turnSpecularOn();
        if (pointlight2.isOn())
            pointlight2.turnSpecularOn();
        if (pointlight3.isOn())
            pointlight3.turnSpecularOn();
        //pointlight4.turnSpecularOn();
        //diffuseToggle = !diffuseToggle;
        //}
    }
    if (glfwGetKey(window, GLFW_KEY_6) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        /*cout << "1 " << pointlight1.isOn() << endl;
        cout << pointlight2.isOn() << endl;
        cout << pointlight3.isOn() << endl;*/
        if (pointlight1.isOn())
            pointlight1.turnSpecularOff();
        if (pointlight2.isOn())
            pointlight2.turnSpecularOff();
        if (pointlight3.isOn())
            pointlight3.turnSpecularOff();
        //pointlight4.turnSpecularOff();
        //diffuseToggle = !diffuseToggle;
        //}
    }
    if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS && glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
    {
        static bool isFullscreen = false;
        static int windowedWidth, windowedHeight, windowedXPos, windowedYPos;

        if (!isFullscreen) {
            glfwGetWindowSize(window, &windowedWidth, &windowedHeight);
            glfwGetWindowPos(window, &windowedXPos, &windowedYPos);
            const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
            glfwSetWindowMonitor(window, glfwGetPrimaryMonitor(), 0, 0, mode->width, mode->height, mode->refreshRate);
        }
        else {
            glfwSetWindowMonitor(window, nullptr, windowedXPos, windowedYPos, windowedWidth, windowedHeight, 0);
        }

        isFullscreen = !isFullscreen;
    }

    

    

}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}


// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}



// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(static_cast<float>(yoffset));
}



unsigned int loadTexture(char const* path, GLenum textureWrappingModeS, GLenum textureWrappingModeT, GLenum textureFilteringModeMin, GLenum textureFilteringModeMax)
{
    unsigned int textureID;
    glGenTextures(1, &textureID);

    int width, height, nrComponents;
    stbi_set_flip_vertically_on_load(true);
    unsigned char* data = stbi_load(path, &width, &height, &nrComponents, 0);
    if (data)
    {
        GLenum format;
        if (nrComponents == 1)
            format = GL_RED;
        else if (nrComponents == 3)
            format = GL_RGB;
        else if (nrComponents == 4)
            format = GL_RGBA;

        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, textureWrappingModeS);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, textureWrappingModeT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, textureFilteringModeMin);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, textureFilteringModeMax);

        stbi_image_free(data);
    }
    else
    {
        std::cout << "Texture failed to load at path: " << path << std::endl;
        stbi_image_free(data);
    }

    return textureID;
}