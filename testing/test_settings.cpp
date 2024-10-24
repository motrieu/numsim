#include <gtest/gtest.h>
#include <src/settings/settings.h>

TEST(SettingsTest, RemoteStartWhiteSpace)
{
    std::string line = "  param = value  # a comment";

    Settings settings;
    settings::removeStartWhitespace(line);

    EXPECT_EQ(linem "param = value  # a comment")
}