#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <variant>
#include <optional>

#include "GCodeCommand.h"

/**
 * @brief The GCodeParser class is used to parse G-code files.
 */
class GCodeParser
{
public:
    /**
     * @brief Parses a file and returns a vector of GCodeCommand objects.
     * @param filename The name of the file to parse.
     * @return A vector of GCodeCommand objects representing the parsed G-code commands.
     */
    std::vector<GCodeCommand> parseFile(const std::string &filename)
    {
        std::ifstream file(filename);
        std::vector<GCodeCommand> commands;

        if (file.is_open())
        {
            std::string line;
            while (std::getline(file, line))
            {
                removeComments(line); // Remove comments from the line
                if (!line.empty())
                {
                    GCodeCommand command = parseLine(line);
                    commands.push_back(command);
                }
            }
            file.close();
        }
        else
        {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }

        return commands;
    }

private:
    /**
     * @brief Removes comments from a line of G-code.
     * @param line The line of G-code to remove comments from.
     */
    void removeComments(std::string &line)
    {
        size_t pos = line.find(';');
        if (pos != std::string::npos)
        {
            line.erase(pos);
        }
    }

    /**
     * @brief Parses a line of G-code and returns a GCodeCommand object.
     * @param line The line of G-code to parse.
     * @return A GCodeCommand object representing the parsed G-code command.
     */
    GCodeCommand parseLine(const std::string &line)
    {
        GCodeCommand command;
        std::istringstream iss(line);
        std::string token;
        iss >> token; // Extract the G-code command (e.g., G1, M104)

        if (!token.empty())
        {
            command.setCode(token);

            while (iss >> token)
            {
                char letter = token[0];
                std::string value = token.substr(1);
                command.addParameter(letter, value);
            }
        }

        return command;
    }
};