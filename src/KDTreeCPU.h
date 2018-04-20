#ifndef KD_TREE_CPU_H
#define KD_TREE_CPU_H

#include <limits>

#include "boundingBox.h"
#include "utils.h"
#include "KDTreeStructs.h"


////////////////////////////////////////////////////
// Constants.
////////////////////////////////////////////////////

const int NUM_TRIS_PER_NODE = 20;
const bool USE_TIGHT_FITTING_BOUNDING_BOXES = false;
//const float INFINITY = std::numeric_limits<float>::max();


////////////////////////////////////////////////////
// KDTreeCPU.
////////////////////////////////////////////////////
class KDTreeCPU
{
public:
	KDTreeCPU( int num_tris, glm::vec3 *tris, int num_verts, glm::vec3 *verts );
	~KDTreeCPU( void );

	// Public traversal method that begins recursive search.
	bool intersect( const glm::vec3 &ray_o, const glm::vec3 &ray_dir, float &t, glm::vec3 &hit_point, glm::vec3 &normal ) const;
	bool intersectNew(const glm::vec3 &ray_o, const glm::vec3 &ray_dir, float &t, glm::vec3 &hit_point, glm::vec3 &normal, int &idx) const;
	bool singleRayStacklessIntersect( const glm::vec3 &ray_o, const glm::vec3 &ray_dir, float &t, glm::vec3 &hit_point, glm::vec3 &normal ) const;

	void buildRopeStructure( void );

	// kd-tree getters.
	KDTreeNode* getRootNode( void ) const;
	int getNumLevels( void ) const;
	int getNumLeaves( void ) const;
	int getNumNodes( void ) const;

	// Input mesh getters.
	int getMeshNumVerts( void ) const;
	int getMeshNumTris( void ) const;
	glm::vec3* getMeshVerts( void ) const;
	glm::vec3* getMeshTris( void ) const;

	// Debug methods.
	void printNumTrianglesInEachNode( KDTreeNode *curr_node, int curr_depth=1 );
	void printNodeIdsAndBounds( KDTreeNode *curr_node );

private:
	// kd-tree variables.
	KDTreeNode *root;
	int num_levels, num_leaves, num_nodes;

	// Input mesh variables.
	int num_verts, num_tris;
	glm::vec3 *verts, *tris;

	KDTreeNode* constructTreeMedianSpaceSplit( int num_tris, int *tri_indices, boundingBox bounds, int curr_depth );

	// Private recursive traversal method.
	bool intersect( KDTreeNode *curr_node, const glm::vec3 &ray_o, const glm::vec3 &ray_dir, float &t, glm::vec3 &normal ) const;
	bool intersectNewRecur(KDTreeNode *curr_node, const glm::vec3 &ray_o, const glm::vec3 &ray_dir, float &t, glm::vec3 &normal, int &idx) const;
	bool singleRayStacklessIntersect( KDTreeNode *curr_node, const glm::vec3 &ray_o, const glm::vec3 &ray_dir, float &t_entry, float &t_exit, glm::vec3 &normal ) const;

	// Rope construction.
	void buildRopeStructure( KDTreeNode *curr_node, KDTreeNode *ropes[], bool is_single_ray_case=false );
	void optimizeRopes( KDTreeNode *ropes[], boundingBox bbox );

	// Bounding box getters.
	SplitAxis getLongestBoundingBoxSide( glm::vec3 min, glm::vec3 max );
	boundingBox computeTightFittingBoundingBox( int num_verts, glm::vec3 *verts );
	boundingBox computeTightFittingBoundingBox( int num_tris, int *tri_indices );

	// Triangle getters.
	float getMinTriValue( int tri_index, SplitAxis axis );
	float getMaxTriValue( int tri_index, SplitAxis axis );
};

#endif