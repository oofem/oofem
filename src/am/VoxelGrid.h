#pragma once

#include <array>
#include <vector>
#include <tuple>
#include <unordered_map>

#include "CSG.h"
#include "field.h"

/**
 * @brief A struct to represent a 3D model.
 */
struct VoxelNode {
    int id               = -1; /**< The OOFEM id of the node. */
    double coords[3]     = { 0.0, 0.0, 0.0 }; /**< The coordinates of the node. */
    double timeActivated = 0.0; /**< The time the node was activated. */
};

/**
 * @brief A struct to represent a single voxel.
 */
struct Voxel {
    int id = -1; /**< The OOFEM id of the voxel. */
    std::array<int, 8> nodes; /**< The OOFEM ids of the nodes of the voxel. */
    double timeActivated = 0.0; /**< The time the voxel was activated. */
    std::vector<std::tuple<double, double> > vofHistory; /**< The time of each vof increment and a corresponding vof value. */

    /**
     * @brief Get the time the voxel was activated.
     * @return The time the voxel was activated.
     */
    double time_activated()
    {
        return this->timeActivated;
    }

    /**
     * @brief Get the volume of fluid (infill percentage) of the voxel.
     * @return The volume of fluid (infill percentage) of the voxel.
     */
    double vof()
    {
        if ( vofHistory.size() == 0 )
            return 0.0;

        return std::get<1>( vofHistory.back() );
    }
    double getVofAt(double time) {
        if (vofHistory.size() == 0)
            return 0.0;
        for (int i = vofHistory.size()-1; i>=0; i--) {
            if (std::get<0>(vofHistory[i]) <= time) {
                return std::get<1>(vofHistory[i]);
            }
        }
        return 0.0;
    }
};

/**
 * @brief A class to represent a 3D grid of voxels.
 */
class VoxelGrid
{
public:
    /**
     * @brief Constructor for VoxelGrid.
     * @param steps The step sizes in each dimension.
     * @param sizes The sizes of the grid in each dimension.
     */
    VoxelGrid( std::array<double, 3> steps, std::array<int, 3> sizes ) :
        m_steps( steps ), m_sizes( sizes )
    {
    }

    /**
     * @brief Get the number of steps in each dimension.
     * @return The number of steps in each dimension.
     */
    std::array<int, 3> sizes()
    {
        return m_sizes;
    }

    /**
     * @brief Get the step sizes in each dimension.
     * @return The step sizes in each dimension.
     */
    std::array<double, 3> steps()
    {
        return m_steps;
    }

    /**
     * @brief Get the 1D index of a voxel in the grid.
     * @param i The x-coordinate of the voxel.
     * @param j The y-coordinate of the voxel.
     * @param k The z-coordinate of the voxel.
     * @return The 1D index of the voxel.
     */
    int get_index( int i, int j, int k )
    {
        return i + m_sizes[0] * ( j + m_sizes[1] * k );
    }

    /**
     * @brief Get the 3D indices of a voxel from its 1D index.
     * @param i The 1D index of the voxel.
     * @return A tuple containing the x, y, and z indices of the voxel.
     */
    std::tuple<int, int, int> get_indices( int i )
    {
        int x = i % m_sizes[0];
        int y = ( i / m_sizes[0] ) % m_sizes[1];
        int z = i / ( m_sizes[0] * m_sizes[1] );

        return std::make_tuple( x, y, z );
    }

    /**
     * @brief Get the 3D indices of a voxel from a point in space.
     * @param point An array containing the x, y, and z coordinates of the point.
     * @return A tuple containing the x, y, and z indices of the voxel.
     */
    std::tuple<int, int, int> get_indices_from_point( std::array<double, 3> point )
    {
        int x = point[0] / m_steps[0];
        int y = point[1] / m_steps[1];
        int z = point[2] / m_steps[2];

        return std::make_tuple( x, y, z );
    }

    /**
     * @brief Get the 3D vector of a voxel (origin) from its 1D index.
     * @param i The 1D index of the voxel.
     * @return An array containing the x, y, and z coordinates of the voxel.
     */
    std::array<double, 3> get_vector_3d( int i )
    {
        auto [x, y, z] = get_indices( i );
        return { x * m_steps[0], y * m_steps[1], z * m_steps[2] };
    }

    std::array<int, 8> giveVoxelNodeIds( int index )
    {
        auto [x, y, z] = get_indices( index );

        std::array<int, 8> nodeIds;
        // top surface
        nodeIds[0] = get_index( x, y, z + 1 );
        nodeIds[1] = get_index( x, y + 1, z + 1 );
        nodeIds[2] = get_index( x + 1, y + 1, z + 1 );
        nodeIds[3] = get_index( x + 1, y, z + 1 );
        // bottom surface
        nodeIds[4] = get_index( x, y, z );
        nodeIds[5] = get_index( x, y + 1, z );
        nodeIds[6] = get_index( x + 1, y + 1, z );
        nodeIds[7] = get_index( x + 1, y, z );
        return nodeIds;
    }

    void activateNodes( int index, double timeActivated )
    {
        auto voxelNodes = giveVoxelNodeIds( index );
        for ( int i = 0; i < 8; i++ ) {
            if ( !is_active_node( voxelNodes[i] ) ) {
                auto coords                          = get_vector_3d( voxelNodes[i] );
                m_nodes[voxelNodes[i]]               = VoxelNode();
                m_nodes[voxelNodes[i]].id            = m_nodes.size();
                m_nodes[voxelNodes[i]].coords[0]     = coords[0];
                m_nodes[voxelNodes[i]].coords[1]     = coords[1];
                m_nodes[voxelNodes[i]].coords[2]     = coords[2];
                m_nodes[voxelNodes[i]].timeActivated = timeActivated;
            }
        }
    }

    /**
     * @brief Create a Vovel at given index and at given time.
     * @param index The index of the voxel to activate.
     * @param timeActivated The time the voxel was activated.
     * 
     */
    Voxel &createVoxel( int index, double timeActivated )
    {
        if ( !is_active( index ) ) {
            // Auto-activate and number nodes
            activateNodes( index, timeActivated );

            // Create new voxel
            m_voxels[index] = Voxel();

            // Calculate one-based OOFEM element id
            m_voxels[index].id = m_voxels.size();

            // Store node ids (OOFEM ids)
            auto ids = giveVoxelNodeIds( index );
            for ( int i = 0; i < 8; i++ ) {
                m_voxels[index].nodes[i] = m_nodes[ids[i]].id;
            }
            // Initialize vof history
            m_voxels[index].vofHistory.push_back( std::make_tuple( 0, 0.0 ) );
            m_voxels[index].timeActivated = timeActivated;
        }

        return m_voxels[index];
    }
    
    /**
     * @brief Activate a voxel at a given index in a given time and store the corresponding vof value.
     * @param index The index of the voxel.
     * @param timeActivated The time the voxel was activated.
     * @param vofIncrement The volume of fluid increment (new material volume [0-1]) in the voxel.
     */
    Voxel &activate( int index, double timeActivated, double vofIncrement )
    {
        this->createVoxel( index, timeActivated );

        double vof_before = m_voxels[index].vof();
        double vof_after  = vof_before + vofIncrement;
        if ( vof_after > 1.0 ) {
            // redistribute the excess material to neighboring voxels
            double excess = vof_after - 1.0;
            vof_after     = 1.0;
            m_voxels[index].vofHistory.push_back( std::make_tuple( timeActivated, vof_after ) );
            // check neighbors in xy plane
            auto [x, y, z] = get_indices( index );
            std::vector<int> neighbor_indices, neighbor_indices_workingset;
            if ( x > 0 )
                neighbor_indices.push_back( get_index( x - 1, y, z ) );
            if ( x < m_sizes[0] - 1 )
                neighbor_indices.push_back( get_index( x + 1, y, z ) );
            if ( y > 0 )        
                neighbor_indices.push_back( get_index( x, y - 1, z ) );
            if ( y < m_sizes[1] - 1 )
                neighbor_indices.push_back( get_index( x, y + 1, z ) );
            if ( x > 0 && y > 0 )
                neighbor_indices.push_back( get_index( x - 1, y - 1, z ) );
            if ( x > 0 && y < m_sizes[1] - 1 )
                neighbor_indices.push_back( get_index( x - 1, y + 1, z ) );
            if ( x < m_sizes[0] - 1 && y > 0 )
                neighbor_indices.push_back( get_index( x + 1, y - 1, z ) );
            if ( x < m_sizes[0] - 1 && y < m_sizes[1] - 1 )
                neighbor_indices.push_back( get_index( x + 1, y + 1, z ) );
            
            double excess_per_neighbor = excess / neighbor_indices.size();
            std::map<int, double> neighborVof; // to accumulate vof increments per neighbor
            for ( auto ni : neighbor_indices ) {
                createVoxel( ni, timeActivated ); // ensure neighbor exists
                neighborVof[ni] = m_voxels[ni].vof();
            }
            neighbor_indices_workingset = neighbor_indices;
            while ( (excess > 0) && (neighbor_indices_workingset.size() > 0 )) {
                excess = 0.0;
                std::vector<int> neighbor_indices2;
                for ( auto ni : neighbor_indices_workingset ) {
                    if (neighborVof[ni] < 1.0) {
                        double vof_after = neighborVof[ni] + excess_per_neighbor;
                        if (vof_after > 1.0) {
                            excess += vof_after - 1.0;
                            vof_after = 1.0;
                            neighborVof[ni] += excess_per_neighbor;
                        } else {
                            neighborVof[ni] = vof_after;
                            neighbor_indices2.push_back( ni ); // remember voxel with free capacity
                        }
                    } else {
                        excess += excess_per_neighbor;
                    }
                }
                excess_per_neighbor = excess / neighbor_indices2.size();
                neighbor_indices_workingset = neighbor_indices2;
            }
            // update neighbors
            for ( auto ni : neighbor_indices ) {
                if (neighborVof[ni] > m_voxels[ni].vof()) {
                    m_voxels[ni].vofHistory.push_back( std::make_tuple( timeActivated, neighborVof[ni] ) );
                }
            }
        } else {
            // vof increment fully accomodated
            m_voxels[index].vofHistory.push_back( std::make_tuple( timeActivated, vof_after ) );
        }   

        return m_voxels[index];
    }

    /**
     * @brief Check if a voxel at a given index is active.
     * @param index The index of the voxel.
     * @return True if the voxel is active, false otherwise.
     */
    bool is_active( int index )
    {
        return m_voxels.find( index ) != m_voxels.end();
    }

    /**
     * @brief Check if a node at a given index is active.
     * @param index The index of the node.
     * @return True if the node is active, false otherwise.
     */
    bool is_active_node( int index )
    {
        return m_nodes.find( index ) != m_nodes.end();
    }

    /**
     * @brief Get the number of active elements (voxels) in the grid.
     * @return The number of active elements.
     */
    int active_elements()
    {
        return m_voxels.size();
    }

    /**
     * @brief Get the number of active nodes in the grid.
     * @return The number of active nodes.
     */
    int active_nodes()
    {
        return m_nodes.size();
    }

    /**
     * @brief Get the size of the grid.
     * @return The size of the grid.
     */
    int size()
    {
        return m_sizes[0] * m_sizes[1] * m_sizes[2];
    }

    /**
     * @brief Get the voxel at a given index.
     * @param index The index of the voxel.
     * @return The voxel at the given index.
     */
    Voxel &get_voxel( int index )
    {
        return m_voxels[index];
    }

    /**
     * @brief Get the node at a given index.
     * @param index The index of the node.
     * @return The node at the given index.
     */
    VoxelNode &get_node( int index )
    {
        return m_nodes[index];
    }

    /**
     * @brief Get the map of nodes in the grid.
     * @return The map of nodes in the grid.
     */
    std::unordered_map<int, Voxel> &giveVoxels()
    {
        return m_voxels;
    }

    /**
     * @brief Get the map of nodes in the grid.
     * @return The map of nodes in the grid.
     */
    std::unordered_map<int, VoxelNode> &giveNodes()
    {
        return m_nodes;
    }

    /**
     * @brief Get the triangulated model of a voxel at a given index.
     * @param index The index of the voxel.
     * @return The model of the voxel at the given index.
     */
    Model get_model( int index )
    {
        auto orig = get_vector_3d( index );
        double dx = m_steps[0];
        double dy = m_steps[1];
        double dz = m_steps[2];

        // std::cout << "C" << orig[0] << " " << orig[1] << " " << orig[2] << std::endl;
        // std::cout << "S" << dx << " " << dy << " " << dz << std::endl;

        Vertex v1( { { orig[0], orig[1], orig[2] }, { -1, -1, -1 } } );
        Vertex v2( { { orig[0] + dx, orig[1], orig[2] }, { 1, -1, -1 } } );
        Vertex v3( { { orig[0] + dx, orig[1] + dy, orig[2] }, { 1, 1, -1 } } );
        Vertex v4( { { orig[0], orig[1] + dy, orig[2] }, { -1, 1, -1 } } );

        Vertex v5( { { orig[0], orig[1], orig[2] + dz }, { -1, -1, 1 } } );
        Vertex v6( { { orig[0] + dx, orig[1], orig[2] + dz }, { 1, -1, 1 } } );
        Vertex v7( { { orig[0] + dx, orig[1] + dy, orig[2] + dz }, { 1, 1, 1 } } );
        Vertex v8( { { orig[0], orig[1] + dy, orig[2] + dz }, { -1, 1, 1 } } );

        Polygon p1( { v1, v2, v6, v5 } ); // bottom
        Polygon p2( { v8, v7, v3, v4 } ); // top
        Polygon p3( { v1, v5, v8, v4 } ); // left
        Polygon p4( { v6, v2, v3, v7 } ); // right
        Polygon p5( { v8, v5, v6, v7 } ); // front
        Polygon p6( { v4, v3, v2, v1 } ); // back

        return modelfrompolygons( { p1, p2, p3, p4, p5, p6 } );
    }

    void writeVTK( std::string filename )
    {
        std::ofstream file( filename );
        file << "# vtk DataFile Version 2.0\n";
        file << "Voxels\n";
        file << "ASCII\n";
        file << "DATASET UNSTRUCTURED_GRID\n";
        file << "POINTS " << m_voxels.size() * 8 << " double\n";

        for ( const auto &[i, voxel] : m_voxels ) {
            auto pt = get_vector_3d( i );

            // Create the 8 points of a cube
            file << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
            file << pt[0] + m_steps[0] << " " << pt[1] << " " << pt[2] << "\n";
            file << pt[0] + m_steps[0] << " " << pt[1] + m_steps[1] << " " << pt[2] << "\n";
            file << pt[0] << " " << pt[1] + m_steps[1] << " " << pt[2] << "\n";

            file << pt[0] << " " << pt[1] << " " << pt[2] + m_steps[2] << "\n";
            file << pt[0] + m_steps[0] << " " << pt[1] << " " << pt[2] + m_steps[2] << "\n";
            file << pt[0] + m_steps[0] << " " << pt[1] + m_steps[1] << " " << pt[2] + m_steps[2] << "\n";
            file << pt[0] << " " << pt[1] + m_steps[1] << " " << pt[2] + m_steps[2] << "\n";
        }

        file << "CELLS " << m_voxels.size() << " " << m_voxels.size() * 9 << "\n";

        for ( size_t i = 0; i < m_voxels.size(); i++ ) {
            file << "8 ";
            for ( int j = 0; j < 8; j++ ) {
                file << i * 8 + j << " ";
            }
            file << "\n";
        }

        file << "CELL_TYPES " << m_voxels.size() << "\n";

        for ( size_t i = 0; i < m_voxels.size(); i++ ) {
            file << "12\n";
        }

        file << "CELL_DATA " << m_voxels.size() << "\n";
        file << "SCALARS vof double 1\n";
        file << "LOOKUP_TABLE default\n";

        for ( auto &[index, voxel] : m_voxels ) {
            file << voxel.vof() << "\n";
        }
    }

private:
    std::unordered_map<int, Voxel> m_voxels; /**< The map of voxels in the grid. */
    std::unordered_map<int, VoxelNode> m_nodes; /**< The map of nodes in the grid. */

    std::array<double, 3> m_steps; /**< The step sizes in each dimension. */
    std::array<int, 3> m_sizes; /**< The sizes of the grid in each dimension. */
};

