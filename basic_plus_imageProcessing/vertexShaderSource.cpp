const char* vertexShaderSource = R"glsl(
#version 330 core
layout (location = 0) in vec3 vertexPosition; // The position variable has attribute position 0
layout (location = 1) in vec3 vertexColor; // The color variable has attribute position 1

uniform mat4 view;
uniform mat4 projection;

out vec3 fragmentColor; // Specify a color output to the fragment shader

void main() {
    gl_Position = projection * view * vec4(vertexPosition, 1.0);
    fragmentColor = vertexColor; // Pass the color to the fragment shader
}
)glsl";
