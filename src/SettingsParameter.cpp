#include "SettingsParameter.h"
#include "Log.h"
#include <iostream>
#include <string>
#include <cstdlib>


SettingsParameter::SettingsParameter(const std::string& name,
    const std::string& defaultValue, UserDefinedStatus userDefined,
    DeprecationStatus deprecated) :
        _name(name), _value(defaultValue), _userDefined(userDefined),
        _isUserDefined(false), _deprecated(deprecated)
{
}


SettingsParameter::SettingsParameter(const SettingsParameter& other) :
    _name(other._name), _value(other._value), _userDefined(other._userDefined),
    _isUserDefined(other._isUserDefined), _deprecated(other._deprecated)
{
}


SettingsParameter& SettingsParameter::operator=(const SettingsParameter &other)
{
    _name = other._name;
    _value = other._value;
    _userDefined = other._userDefined;
    _isUserDefined = other._isUserDefined;
    _deprecated = other._deprecated;

    return *this;
}


void SettingsParameter::exitWithErrorWrongType() const
{
    exitWithError("ERROR: Parameter " + _name + " has the wrong type.\n" +
                 "Fix by assigning the parameter a value of the right type.");
}
