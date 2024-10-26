#include <gtest/gtest.h>
#include "settings/settings.h"

TEST(SettingsTest, RemoveStartWhiteSpace)
{
    std::string line = "  param = value  # a comment";

    Settings settings;
    settings.removeStartWhitespace(line);

    EXPECT_EQ(line, "param = value  # a comment");
}