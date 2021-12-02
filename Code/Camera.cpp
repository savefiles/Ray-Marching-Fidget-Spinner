#include "Camera.h"

Camera::Camera(HWND* window) : Camera(glm::vec3(0.0f, 3.0f, -3.0f), window) {}

Camera::Camera(glm::vec3 position, HWND* window)
{
	this->window = window;

	// Update the screen center using the HWND.
	UpdateScreenCenter();

	// Set the default values.
	this->position = position;

	forward = glm::normalize(-position);
	up = glm::vec3(0.0f, 1.0f, 0.0f);
	right = glm::normalize(glm::cross(forward, up));
}

void Camera::MoveLocal(glm::vec3 localMoveVector)
{
	// We apply each component of the local move vector individually.
	position += right * localMoveVector.x;
	position += glm::vec3(0.0f, localMoveVector.y, 0.0f);
	position += (forward * glm::vec3(1.0f, 0.0f, 1.0f)) * localMoveVector.z ;
}

void Camera::SetPosition(glm::vec3 worldPoint)
{
	position = worldPoint;
}

void Camera::UpdateMouseOrientation()
{
	// Get the mouse position in pixel coords.
	POINT pt;
	GetCursorPos(&pt);
	ScreenToClient(*window, &pt);

	// Get the mouse delta (change in angle) for both directions.
	float dx = (pt.x - centerOfScreen.x) * mouseSens;
	float dy = (pt.y - centerOfScreen.y) * mouseSens;

	// Limit the pitch of the camera.
	float dot = glm::dot(forward, glm::vec3(0.0f, 1.0f, 0.0f));
	if (dot >=  0.95f && dy <= 0.0f) { 
		dy = 0.0f; 
	}
	if (dot <= -0.95f && dy >= 0.0f) { 
		dy = 0.0f; 
	}

	//if (forward.y >= 0.95f && dy >= 0.0f) { dy = 0; }
	//if (forward.y <= -0.95f && dy <= 0.0f) { dy = 0; }

	// Create and apply rotation quaternion.
	glm::quat rotation = glm::angleAxis(dy, right) * glm::angleAxis(dx, up);
	forward = forward * rotation;

	// Normalize all three camera vectors.
	forward = glm::normalize(forward);
	right = glm::normalize(glm::cross(forward, up));
	// (dont need to normalize up because it doesn't change)

	// Set the cursor position to the center of the screen.
	ClientToScreen(*window, &centerOfScreen);
	SetCursorPos(centerOfScreen.x, centerOfScreen.y);
	ScreenToClient(*window, &centerOfScreen);
}

void Camera::UpdateScreenCenter()
{
	// Get the center of the screen in client coordinates (relative to top left of client area).
	RECT dim;
	GetClientRect(*window, &dim);
	centerOfScreen = { (dim.left + dim.right) / 2, (dim.top + dim.bottom) / 2 };
}

glm::mat4 Camera::GetViewMatrix()
{
	return glm::lookAt(position, position + forward, up);
}

const glm::vec3& Camera::GetPosition()
{
	return position;
}

const glm::vec3& Camera::GetForward()
{
	return forward;
}
