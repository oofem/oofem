#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <variant>
#include <optional>

/**
 * @brief The GCodeCommand class represents a single G-code command.
 */
class GCodeCommand
{
public:
    /**
     * @brief Default constructor for the GCodeCommand class.
     */
    GCodeCommand() = default;

    /**
     * @brief Constructor for the GCodeCommand class.
     * @param cmd The G-code command.
     */
    GCodeCommand(const std::string &cmd) : code(cmd) {}

    /**
     * @brief Method to add a parameter to the G-code command.
     * @param letter The parameter letter.
     * @param value The parameter value.
     */
    void addParameter(char letter, const std::string &value)
    {
        parameters[letter] = value;
    }

    /**
     * @brief Set the G-code command.
     * @param code The G-code command.
     */
    void setCode(const std::string code) { this->code = code; }

    /**
     * @brief Get the G-code command.
     * @return The G-code command.
     */
    std::string getCode() const { return code; }

    /**
     * @brief Get the parameters associated with the G-code command.
     * @return The parameters as a map of parameter letter to parameter value.
     */
    std::map<char, std::string> getParameters() const { return parameters; }

    /**
     * @brief Get the value of a parameter as a string.
     * @param paramLetter The parameter letter.
     * @return The parameter value as a string, or std::nullopt if the parameter is not found.
     */
    std::optional<std::string> param_string(char paramLetter) const
    {
        auto it = parameters.find(paramLetter);
        if (it != parameters.end())
        {
            return it->second;
        }
        return std::nullopt;
    }

    /**
     * @brief Get the value of a parameter as an integer.
     * @param paramLetter The parameter letter.
     * @return The parameter value as an integer, or std::nullopt if the parameter is not found or the conversion fails.
     */
    std::optional<int> param_int(char paramLetter) const
    {
        auto strValue = param_string(paramLetter);
        if (strValue)
        {
            try
            {
                return std::stoi(*strValue);
            }
            catch (...)
            {
                return std::nullopt; // Conversion failed
            }
        }
        return std::nullopt; // Parameter not found
    }

    /**
     * @brief Get the value of a parameter as a double.
     * @param paramLetter The parameter letter.
     * @return The parameter value as a double, or std::nullopt if the parameter is not found or the conversion fails.
     */
    std::optional<double> param_double(char paramLetter) const
    {
        auto strValue = param_string(paramLetter);
        if (strValue)
        {
            try
            {
                return std::stod(*strValue);
            }
            catch (...)
            {
                return std::nullopt; // Conversion failed
            }
        }
        return std::nullopt; // Parameter not found
    }

private:
    std::string code;                       // The G-code command
    std::map<char, std::string> parameters; // The parameters associated with the G-code command
};