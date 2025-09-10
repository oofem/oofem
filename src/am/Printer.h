#pragma once

#include <array>
#include <iostream>
#include <map>
#include <functional>
#include <cmath>
#include <queue>

#include "GCodeCommand.h"
#include "CSG.h"
#include "VoxelGrid.h"

/**
 * @brief The PrintStatistics struct is used to track statistics of the printer.
 */
struct PrintStatistics {
    double distance_moved    = 0; /**< The total distance moved by the printer. */
    double filament_extruded = 0; /**< The total amount of filament extruded by the printer. */
    double time              = 0; /**< The total time taken by the printer to process commands. */
    int    nunber_of_commands = 0; /**< The total number of commands processed by the printer. */
};

/**
 * @brief The ExtruderPositioning enum class represents the positioning mode of the extruder.
 */
enum class ExtruderPositioning {
    Absolute, /**< The extruder is positioned absolutely. */
    Relative /**< The extruder is positioned relatively. */
};

/**
 * @brief The VelocityModel enum class represents the velocity model of the printer.
 * Acceleration and deceleration are handled differently for each model.
 */
enum class VelocityModel {
    Constant, /**< Constant velocity model (infinite acceleration). */
    Trapezoidal, /**< Trapezoidal velocity model (constant acceleration). */
    SCurve /**< S-curve velocity model. */
};

/**
 * @brief The LayerHeight enum class represents the layer height of the printer.
 * The layer height can be constant, first then constant, or automatic.
 */
enum class LayerHeightModel {
    Constant,
    FirstThenConstant,
    Automatic
};

/**
 * @brief The PrinterOptions struct is used to configure the Printer class.
 */
struct PrinterOptions {
    VelocityModel velocityModel = VelocityModel::Trapezoidal; /**< The velocity model of the printer. */

    std::array<double, 3> steps = { 1.0, 1.0, 1.0 }; /**< The step sizes in each dimension. */
    std::array<int, 3> sizes    = { 1000, 1000, 1000 }; /**< The sizes of the grid in each dimension. */

    LayerHeightModel layerHeightModel = LayerHeightModel::Constant; /**< The layer height model of the printer. */
    double layerHeight                = 0.2; /**< The layer height value. */
    double extrusionWidth             = 0.4 * 1.2;
    double depositionTemperature      = 235.; /** Temperature of deposited material */
    double heatBedTemperature         = 60.; /* heat bed temperature */
    double chamberTemperature         = 30.; /* air temperature in printer chamber*/
    double heatTransferFilmCoefficient= 10.;  /* film coefficient for heat transfer between deposited material and air*/
    double depositedMaterialHeatPower = 4200000.0; /* power = specificHeat * density */

    double filamentDiameter = 1.75; /**< The diameter of the filament. Not used*/
};

/**
 * @brief The Printer class is used to process G-code commands and manage the printer state.
 */
class Printer
{
public:
    // Define a type for the callback function
    using CommandCallback = std::function<void( const GCodeCommand & )>;

    /**
     * @brief Constructor for Printer using PrinterOptions.
     */
    Printer( PrinterOptions options )
    {
        velocity_model     = options.velocityModel;
        voxelGrid          = { { options.steps[0], options.steps[1], options.steps[2] }, { options.sizes[0], options.sizes[1], options.sizes[2] } };
        layer_height_model = options.layerHeightModel;
        layer_height       = options.layerHeight;
        extrusion_width    = options.extrusionWidth;

        filament_diameter  = options.filamentDiameter;

        depositionTemperature = options.depositionTemperature; 
        heatBedTemperature    = options.heatBedTemperature;
        chamberTemperature    = options.chamberTemperature;
        heatTransferFilmCoefficient= options.heatTransferFilmCoefficient;
        depositedMaterialHeatPower = options.depositedMaterialHeatPower;

        registerCallbacks();
    }

    /**
     * @brief Get the printer statistics.
     * @return The printer statistics.
     */
    PrintStatistics getStatistics() const
    {
        return statistics;
    }

    /**
     * @brief Register a callback function for a specific G-code command.
     * @param commandCode The G-code command code.
     * @param callback The callback function to be registered.
     */
    void registerCallback( const std::string &commandCode, CommandCallback callback )
    {
        commandCallbacks[commandCode] = callback;
    }

    /**
     * @brief Add a G-code command to the command queue.
     * @param command The G-code command to be added.
     */
    void addCommandToQueue( const GCodeCommand &command )
    {
        commandQueue.push_back( command );
    }

    /**
     * @brief Remove the first G-code command from the command queue.
     */
    void popCommandFromQueue()
    {
        if ( !commandQueue.empty() )
            commandQueue.erase( commandQueue.begin() );
    }

    /**
     * @brief Process a G-code command.
     * @param command The G-code command to be processed.
     */
    void processCommand( const GCodeCommand &command )
    {
        auto it = commandCallbacks.find( command.getCode() );
        if ( it != commandCallbacks.end() ) {
            // Execute the callback associated with the command code
            this->statistics.nunber_of_commands++;
            it->second( command );
        } else {
            // Handle unknown command
            // std::cout << "Unknown command: " << command.getCode() << std::endl;
        }
    }

    VoxelGrid &getGrid()
    {
        return voxelGrid;
    }

    /**
     * @brief Calculate the new position of the printer based on the given values.
     * @param position The current position of the printer.
     * @param xValue The optional X-axis value.
     * @param yValue The optional Y-axis value.
     * @param zValue The optional Z-axis value.
     * @return The new position of the printer.
     */
    std::array<double, 4> calculatePosition( std::array<double, 4> position,
        const std::optional<double> &xValue,
        const std::optional<double> &yValue,
        const std::optional<double> &zValue )
    {
        // If relative position, move by delta
        if ( positioning == ExtruderPositioning::Relative ) {
            if ( xValue.has_value() )
                position[0] += *xValue;
            if ( yValue.has_value() )
                position[1] += *yValue;
            if ( zValue.has_value() )
                position[2] += *zValue;
        }
        // If absolute position, move to the position
        else {
            if ( xValue.has_value() )
                position[0] = *xValue;
            if ( yValue.has_value() )
                position[1] = *yValue;
            if ( zValue.has_value() )
                position[2] = *zValue;
        }

        return position;
    }

    /**
     * @brief Calculate the director vector between two positions.
     * @param position_prev The previous position.
     * @param position_next The next position.
     * @return The director vector.
     */
    std::array<double, 3>
    calculateDirectorVector( const std::array<double, 4> &position_prev,
        const std::array<double, 4> &position_next )
    {
        std::array<double, 3> director_vector = { position_next[0] - position_prev[0],
            position_next[1] - position_prev[1],
            position_next[2] - position_prev[2] };

        // Normalize
        double norm = std::sqrt( director_vector[0] * director_vector[0] + director_vector[1] * director_vector[1] + director_vector[2] * director_vector[2] );

        if ( std::abs( norm ) > 1e-16 ) {
            director_vector[0] /= norm;
            director_vector[1] /= norm;
            director_vector[2] /= norm;
        }

        return director_vector;
    }

    /**
     * @brief Calculate the time taken to move a certain distance with given velocities.
     * @param distance The distance to be moved.
     * @param vs The starting velocity.
     * @param v The target velocity.
     * @param ve The ending velocity.
     * @return The time taken to move the distance.
     */
    double calculateMoveTime( double distance, double vs, double v, double ve )
    {
        double time_acc     = ( v - vs ) / max_acceleration[0];
        double distance_acc = time_acc * ( v - vs ) / 60 / 2;

        double time_decc     = ( v - ve ) / max_acceleration[0];
        double distance_decc = time_decc * ( v - ve ) / 60 / 2;

        return time_acc + time_decc + ( distance - distance_acc - distance_decc ) / v;
    }

    const Model get_model( std::array<double, 4> start, std::array<double, 4> end, double h, double w )
    {
        std::array<double, 3> dir = { end[0] - start[0], end[1] - start[1], 0 };
        std::array<double, 3> n   = { -dir[1], dir[0], dir[2] };

        // Normalize dir
        double norm_dir = std::sqrt( dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2] );
        dir[0] /= norm_dir;
        dir[1] /= norm_dir;
        dir[2] /= norm_dir;

        // Normalize n
        double norm = std::sqrt( n[0] * n[0] + n[1] * n[1] + n[2] * n[2] );
        n[0] *= w / 2 / norm;
        n[1] *= w / 2 / norm;
        n[2] *= w / 2 / norm;

        Vertex v1( { { start[0] + n[0], start[1] + n[1], start[2] + n[2] - h }, { -dir[0] + n[0], -dir[1] + n[1], -1 } } );
        Vertex v2( { { end[0] + n[0], end[1] + n[1], end[2] + n[2] - h }, { dir[0] + n[0], dir[1] + n[1], -1 } } );
        Vertex v3( { { end[0] - n[0], end[1] - n[1], end[2] - n[2] - h }, { dir[0] - n[0], dir[1] - n[1], -1 } } );
        Vertex v4( { { start[0] - n[0], start[1] - n[1], start[2] - n[2] - h }, { -dir[0] - n[0], -dir[1] - n[1], -1 } } );

        Vertex v5( { { start[0] + n[0], start[1] + n[1], start[2] + n[2] }, { -dir[0] + n[0], -dir[1] + n[1], 1 } } );
        Vertex v6( { { end[0] + n[0], end[1] + n[1], end[2] + n[2] }, { dir[0] + n[0], dir[1] + n[1], 1 } } );
        Vertex v7( { { end[0] - n[0], end[1] - n[1], end[2] - n[2] }, { dir[0] - n[0], dir[1] - n[1], 1 } } );
        Vertex v8( { { start[0] - n[0], start[1] - n[1], start[2] - n[2] }, { -dir[0] - n[0], -dir[1] - n[1], 1 } } );

        // Polygon p1({v4, v3, v2, v1}); // bottom
        // Polygon p2({v5, v6, v7, v8}); // top
        // Polygon p3({v4, v8, v7, v3}); // left
        // Polygon p4({v1, v2, v6, v5}); // right
        // Polygon p5({v4, v1, v5, v8}); // front
        // Polygon p6({v2, v3, v7, v6}); // back

        Polygon p1( { v1, v2, v3, v4 } ); // bottom
        Polygon p2( { v8, v7, v6, v5 } ); // top
        Polygon p3( { v3, v7, v8, v4 } ); // left
        Polygon p4( { v5, v6, v2, v1 } ); // right
        Polygon p5( { v8, v5, v1, v4 } ); // front
        Polygon p6( { v6, v7, v3, v2 } ); // back

        return modelfrompolygons( { p1, p2, p3, p4, p5, p6 } );
    }

    /**
     * @brief Handle the G1 command to manage printer positions.
     * @param command The G-code command.
     */
    void handleG1Command( const GCodeCommand &command )
    {

        // Retrieve move
        auto xValue = command.param_double( 'X' );
        auto yValue = command.param_double( 'Y' );
        auto zValue = command.param_double( 'Z' );
        auto fValue = command.param_double( 'F' );


        // Handle new desired feedrate
        if ( fValue.has_value() ) {
            feedrate = *fValue;
        }


        // Previous position
        std::array<double, 4> prev_position = { position[0], position[1], position[2], position[3] };
        position                            = calculatePosition( position, xValue, yValue, zValue );


        // Calculate distance moved
        double dx       = position[0] - prev_position[0];
        double dy       = position[1] - prev_position[1];
        double dz       = position[2] - prev_position[2];
        double distance = std::sqrt( dx * dx + dy * dy + dz * dz );
        double dxy      = std::sqrt( dx * dx + dy * dy );
        statistics.distance_moved += distance;

        std::array<double, 4> tmp_position = position;
        double lim_feedrate                = feedrate / 60;


        // Look into command queue to calculate target velocity at the end of the current move
        for ( size_t i = 0; i < commandQueue.size(); i++ ) {
            auto newPosition = calculatePosition( tmp_position, commandQueue[i].param_double( 'X' ), commandQueue.front().param_double( 'Y' ), commandQueue.front().param_double( 'Z' ) );
            auto dirVector   = calculateDirectorVector( tmp_position, newPosition );

            std::array<double, 3> jerk = { dirVector[0] * feedrate / 60,
                dirVector[1] * feedrate / 60,
                dirVector[2] * feedrate / 60 };

            if ( jerk[0] > max_jerk[0] ) {
                lim_feedrate = std::min( lim_feedrate, jerk[0] );
            }

            if ( jerk[1] > max_jerk[1] ) {
                lim_feedrate = std::min( lim_feedrate, jerk[1] );
            }

            if ( jerk[2] > max_jerk[2] ) {
                lim_feedrate = std::min( lim_feedrate, jerk[2] );
            }

            tmp_position = newPosition;

            if ( commandQueue[i].getCode() == "G1" )
                break;
        }


        // Calculate print time of the current move
        double move_time = calculateMoveTime( distance, real_feedrate / 60, feedrate / 60, lim_feedrate );
        statistics.time += move_time;

        // Set real velocity at the end of the current move
        real_feedrate = lim_feedrate * 60;


        // Handle extrusion
        auto eValue = command.param_double( 'E' );
        if ( eValue.has_value() ) {

            // If relative extruder positioning, extrude by delta
            if ( extruder_positioning == ExtruderPositioning::Relative ) {
                position[3] += *eValue;
            } else // If absolute extruder positioning, extrude to the position
            {
                position[3] = *eValue;
            }

            // Calculate filament extruded
            double de = position[3] - prev_position[3];
            statistics.filament_extruded += de;

            // Activate elements if printing
            if ( dxy > 0 && de > 0 && dz == 0 ) {
                // Activate grid elements
                double h   = layer_height_model == LayerHeightModel::Constant ? layer_height : 0.2;
#if 0                
                double r   = filament_diameter / 2;
                double w   = M_PI * r * r * de / ( dxy * h );
#else
                double w = extrusion_width;
#endif

                Model move = get_model( prev_position, position, h, w );

                auto pt  = voxelGrid.get_indices_from_point( { prev_position[0], prev_position[1], prev_position[2] } );
                auto pt2 = voxelGrid.get_indices_from_point( { position[0], position[1], position[2] - h } );
                int top  = std::get<2>( pt ) + 1;
                int bot  = std::get<2>( pt2 ) - 1;

                int minx = std::min( std::get<0>( pt ), std::get<0>( pt2 ) ) - 1;
                int maxx = std::max( std::get<0>( pt ), std::get<0>( pt2 ) ) + 1;
                int miny = std::min( std::get<1>( pt ), std::get<1>( pt2 ) ) - 1;
                int maxy = std::max( std::get<1>( pt ), std::get<1>( pt2 ) ) + 1;

                // Volume of single voxel
                double elvol = getGrid().get_model( 0 ).volume();

                std::vector<int> indicesToCheck;

                for ( int z = bot; z <= top; z++ ) {
                    for ( int x = minx; x <= maxx; x++ ) {
                        for ( int y = miny; y <= maxy; y++ ) {
                            indicesToCheck.push_back( voxelGrid.get_index( x, y, z ) );
                        }
                    }
                }


                std::vector<std::tuple<int, double, double> > intersected_indices;
#ifdef _OPENMP
#pragma omp parallel
#endif
                {
                    std::vector<std::tuple<int, double, double> > intersected_indices_local;
#ifdef _OPENMP
#pragma omp for nowait
#endif
                    for ( size_t ix = 0; ix < indicesToCheck.size(); ix++ ) {
                        int i          = indicesToCheck[ix];
                        Model voxModel = voxelGrid.get_model( i );
                        Model is       = csgintersection( move, voxModel );

                        if ( is.volume() > 0.0 ) {
                            double vofFrac = is.volume() / elvol;
                            intersected_indices_local.push_back( std::make_tuple( i, statistics.time, vofFrac ) );
                        }
                    }
#ifdef _OPENMP
#pragma omp critical
#endif
                    intersected_indices.insert( intersected_indices.end(), intersected_indices_local.begin(), intersected_indices_local.end() );
                }

                for ( const auto &v : intersected_indices ) {
                    voxelGrid.activate( std::get<0>( v ), std::get<1>( v ), std::get<2>( v ) );
                }
            }
        }
        //printf("[%d] G1: Ts=%.2f, Te=%.2f, dist=%.2f, h=%.2f, E=%.2f\n", this->statistics.nunber_of_commands, statistics.time-move_time, statistics.time, distance, prev_position[2], eValue.value_or(0.0));
    }

    /**
     * @brief Handle the M83 command to set the extruder positioning mode to relative.
     * @param command The G-code command.
     */
    void handleM83Command( const GCodeCommand &command )
    {
        extruder_positioning = ExtruderPositioning::Relative;
    }

    /**
     * @brief Constructor to register default callbacks.
     */
    Printer()
    {
        registerCallbacks();
    }

private:
    void registerCallbacks()
    {
        registerCallback( "G1", std::bind( &Printer::handleG1Command, this, std::placeholders::_1 ) );
        registerCallback( "M83", std::bind( &Printer::handleM83Command, this, std::placeholders::_1 ) );
    }

    VoxelGrid voxelGrid = { { 0.5, 0.5, 0.5 }, { 1000, 1000, 1000 } }; /**< 3d Voxel Grid */

    std::map<std::string, CommandCallback> commandCallbacks; /**< Map of command codes to callback functions. */
    std::vector<GCodeCommand> commandQueue; /**< Queue of G-code commands to be processed. */
    PrintStatistics statistics; /**< The printer statistics. */

    // Machine positioning
    std::array<double, 4> position           = { 0.0, 0.0, 0.0, 0.0 }; /**< The position of the print head in X, Y, Z, and E axes. */
    ExtruderPositioning positioning          = ExtruderPositioning::Absolute; /**< The positioning mode of the printhead in space (XYZ). */
    ExtruderPositioning extruder_positioning = ExtruderPositioning::Absolute; /**< The positioning mode of the extruder (E axis). */

    // Machine velocity
    VelocityModel velocity_model = VelocityModel::Trapezoidal; /**< The velocity model of the printer. */
    double feedrate              = 1000.0; /**< The current desired velocity of the printer. */
    double real_feedrate         = 0.0; /**< The real velocity of the printer. */

    // Machine limits
    double max_acceleration[4] = { 500.0, 500.0, 100.0, 5000.0 }; /**< The maximum acceleration for X, Y, Z, and E axes. */
    double max_jerk[4]         = { 8.0, 8.0, 0.4, 5.0 };

    // Dimensions
    LayerHeightModel layer_height_model = LayerHeightModel::Constant; /**< The layer height model of the printer. */
    double layer_height                 = 0.2; /**< The layer height value. */
    double nozzle_diameter              = 0.4; /**< Nozzle diameter, not used */
    double filament_diameter            = 1.75; /**< Filament diameter, not used */
    double extrusion_multiplier         = 1.0; /**< Extrusion multiplier */
    double extrusion_width              = nozzle_diameter * 1.2; /**< Current extrusion width */
public:
    double depositionTemperature      = 235.; /** Temperature of deposited material */
    double heatBedTemperature         = 60.; /* heat bed temperature */
    double chamberTemperature         = 30.; /* air temperature in printer chamber*/
    double heatTransferFilmCoefficient= 10.;  /* film coefficient for heat transfer between deposited material and air*/
    double depositedMaterialHeatPower = 4200000.0; /* power = specificHeat * density */
};