#pragma once

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <windows.h>

class Camera
{
public:

	// Constructors.
	Camera(HWND* window);
	Camera(glm::vec3 position, HWND* window);

	// Move the camera based on local space.
	void MoveLocal(glm::vec3 localMoveVector);

	// Set the position of the camera.
	void SetPosition(glm::vec3 worldPoint);

	// Rotate the camera based on mouse position.
	void UpdateMouseOrientation();

	// Change the center of screen var. Call on window move/resize.
	void UpdateScreenCenter();

	// Getters (ish)
	glm::mat4 GetViewMatrix();
	const glm::vec3& GetPosition();
	const glm::vec3& GetForward();

private:

	HWND* window;

	// Represent the camera.
	glm::vec3 position;
	
	// 3 vectors representing the camera's orientation.
	glm::vec3 forward;
	glm::vec3 up;
	glm::vec3 right;

	// Half screen size.
	POINT centerOfScreen;

	// Mouse rotation sensitivity.
	float mouseSens = 0.003f;
};

