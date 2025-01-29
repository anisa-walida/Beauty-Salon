#ifndef CYLINDER_H
#define CYLINDER_H

#include <glad/glad.h>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "shader.h"

#define PI 3.1416

class Cylinder {
public:
    glm::vec3 ambient;
    glm::vec3 diffuse;
    glm::vec3 specular;
    float shininess;

    // Constructor
    Cylinder(float radius = 1.0f, float height = 2.0f, int sectorCount = 36,
        glm::vec3 amb = glm::vec3(0.5f, 0.5f, 1.0f),
        glm::vec3 diff = glm::vec3(0.8f, 0.8f, 1.0f),
        glm::vec3 spec = glm::vec3(1.0f, 1.0f, 1.0f), float shiny = 32.0f)
        : verticesStride(32) {
        set(radius, height, sectorCount, amb, diff, spec, shiny);
        buildCoordinatesAndIndices();
        buildVertices();

        // Create VAO for textured cylinder
        glGenVertexArrays(1, &cylinderVAO);
        glBindVertexArray(cylinderVAO);

        unsigned int VBO, EBO;
        glGenBuffers(1, &VBO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, getVertexSize(), getVertices(), GL_STATIC_DRAW);

        glGenBuffers(1, &EBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, getIndexSize(), getIndices(), GL_STATIC_DRAW);

        glEnableVertexAttribArray(0);  // Position
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, verticesStride, (void*)0);

        glEnableVertexAttribArray(1);  // Normal
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, verticesStride, (void*)(3 * sizeof(float)));

        glEnableVertexAttribArray(2);  // Texture coordinates
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, verticesStride, (void*)(6 * sizeof(float)));

        glBindVertexArray(0); // Unbind VAO
    }

    ~Cylinder() {}

    void set(float radius, float height, int sectors, glm::vec3 amb, glm::vec3 diff, glm::vec3 spec, float shiny) {
        this->radius = radius;
        this->height = height;
        this->sectorCount = (sectors < 3) ? 3 : sectors;
        this->ambient = amb;
        this->diffuse = diff;
        this->specular = spec;
        this->shininess = shiny;
    }

    unsigned int getVertexSize() const { return (unsigned int)vertices.size() * sizeof(float); }
    const float* getVertices() const { return vertices.data(); }
    unsigned int getIndexSize() const { return (unsigned int)indices.size() * sizeof(unsigned int); }
    const unsigned int* getIndices() const { return indices.data(); }
    unsigned int getIndexCount() const { return (unsigned int)indices.size(); }

    // Draw the textured cylinder
    void drawCylinderTexture(Shader& shader, glm::mat4 model, unsigned int diffuseMap, unsigned int specularMap) const {
        shader.use();
        shader.setInt("material.diffuse", 0);
        shader.setInt("material.specular", 1);
        shader.setFloat("material.shininess", shininess);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, diffuseMap);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, specularMap);

        glBindVertexArray(cylinderVAO);
        shader.setMat4("model", model);
        glDrawElements(GL_TRIANGLES, getIndexCount(), GL_UNSIGNED_INT, (void*)0);
        glBindVertexArray(0);
    }

private:
    void buildCoordinatesAndIndices() {
        float x, y, z, s, t;
        float sectorStep = 2 * PI / sectorCount;

        // Bottom circle center vertex
        coordinates.push_back(0.0f);  // x
        coordinates.push_back(-height / 2);  // y
        coordinates.push_back(0.0f);  // z
        normals.push_back(0.0f);  // nx
        normals.push_back(-1.0f); // ny
        normals.push_back(0.0f);  // nz
        texCoords.push_back(0.5f);  // s
        texCoords.push_back(0.5f);  // t

        // Bottom circle perimeter
        for (int i = 0; i <= sectorCount; ++i) {
            float angle = i * sectorStep;
            x = radius * cosf(angle);
            z = radius * sinf(angle);
            y = -height / 2;

            coordinates.push_back(x);
            coordinates.push_back(y);
            coordinates.push_back(z);

            normals.push_back(0.0f);
            normals.push_back(-1.0f);
            normals.push_back(0.0f);

            s = (cosf(angle) + 1.0f) * 0.5f;
            t = (sinf(angle) + 1.0f) * 0.5f;
            texCoords.push_back(s);
            texCoords.push_back(t);
        }

        // Top circle center vertex
        coordinates.push_back(0.0f);  // x
        coordinates.push_back(height / 2);  // y
        coordinates.push_back(0.0f);  // z
        normals.push_back(0.0f);  // nx
        normals.push_back(1.0f);  // ny
        normals.push_back(0.0f);  // nz
        texCoords.push_back(0.5f);  // s
        texCoords.push_back(0.5f);  // t

        // Top circle perimeter
        for (int i = 0; i <= sectorCount; ++i) {
            float angle = i * sectorStep;
            x = radius * cosf(angle);
            z = radius * sinf(angle);
            y = height / 2;

            coordinates.push_back(x);
            coordinates.push_back(y);
            coordinates.push_back(z);

            normals.push_back(0.0f);
            normals.push_back(1.0f);
            normals.push_back(0.0f);

            s = (cosf(angle) + 1.0f) * 0.5f;
            t = (sinf(angle) + 1.0f) * 0.5f;
            texCoords.push_back(s);
            texCoords.push_back(t);
        }

        // Indices for the bottom circle
        int centerIndexBottom = 0;
        for (int i = 1; i <= sectorCount; ++i) {
            indices.push_back(centerIndexBottom);
            indices.push_back(i);
            indices.push_back(i + 1);
        }

        // Indices for the top circle
        int centerIndexTop = (int)coordinates.size() / 3 - (sectorCount + 2);
        for (int i = 1; i <= sectorCount; ++i) {
            indices.push_back(centerIndexTop);
            indices.push_back(centerIndexTop + i);
            indices.push_back(centerIndexTop + i + 1);
        }

        // Indices for the sides
        for (int i = 1; i <= sectorCount; ++i) {
            int next = i + 1;

            // First triangle
            indices.push_back(i);
            indices.push_back(i + sectorCount + 1);
            indices.push_back(next);

            // Second triangle
            indices.push_back(next);
            indices.push_back(i + sectorCount + 1);
            indices.push_back(next + sectorCount + 1);
        }
    }

    void buildVertices() {
        size_t count = coordinates.size();
        for (size_t i = 0; i < count; i += 3) {
            // Position
            vertices.push_back(coordinates[i]);
            vertices.push_back(coordinates[i + 1]);
            vertices.push_back(coordinates[i + 2]);

            // Normal
            vertices.push_back(normals[i]);
            vertices.push_back(normals[i + 1]);
            vertices.push_back(normals[i + 2]);

            // Texture Coordinates
            vertices.push_back(texCoords[i / 3 * 2]);
            vertices.push_back(texCoords[i / 3 * 2 + 1]);
        }
    }

    unsigned int cylinderVAO;
    float radius;
    float height;
    int sectorCount;
    std::vector<float> vertices;
    std::vector<float> coordinates;
    std::vector<float> normals;
    std::vector<float> texCoords;
    std::vector<unsigned int> indices;
    const int verticesStride;
};

#endif /* CYLINDER_H */