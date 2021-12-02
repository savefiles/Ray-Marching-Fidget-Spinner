#pragma once
#include "Pipe.h"

class PipeRay : public Pipe
{
public:
	PipeRay(VkDevice dev, VkDescriptorSetLayout desc_layout, VkRenderPass render_pass);
	~PipeRay();
};

