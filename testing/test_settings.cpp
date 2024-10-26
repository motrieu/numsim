#include <gtest/gtest.h>
#include "settings/settings.h"

class SettingsTest : public testing::Test
{
protected:
    Settings settings;
    std::string line = "  param = value  # a comment";
};

TEST_F(SettingsTest, RemoveStartWhiteSpace)
{
    settings.removeStartWhitespace(line);

    EXPECT_EQ(line, "param = value  # a comment");
}

TEST_F(SettingsTest, SetParanter)
{
    std::string paramName = "re";
    std::string paramVal = "5";

    settings.setParameter(paramName, paramVal);

    EXPECT_EQ(settings.re, 5);
}

TEST_F(SettingsTest, ExtractValueString)
{
    std::string valueString = settings.extractValueString(line);

    EXPECT_EQ(valueString, "value");
}

TEST_F(SettingsTest, ExtractParameterName)
{
    std::string valueString = settings.extractParameterName(line);
    EXPECT_EQ(valueString, "");

    line = "param = value  # a comment";
    valueString = settings.extractParameterName(line);

    EXPECT_EQ(valueString, "param");
}