using Test, Markdown

using .MDDocTestTools

function test_doctest_tools()
    @testset verbose = true "Section text extraction" begin
        text = """
          # Header 1
          This is some text
          ## Header 2
          This is some other text.

          On more paragraphs.
          This is a line in a paragraph.

          ## Header 2.5
          A list
          - one
          - two
          - three

          ## Header 2.6
          A list
          * one
          * two
          * three

          ## A link

          [this is a link](to somewhere far away)

          ## Some code

          ```C
          int x,y
          ```

          ## Header 3
          This is the end bit.

          Which is quite long and has two paragraphs
          (and is at the end).
          """

        test_md = Markdown.parse(text)

        function test_section_extraction(section_name, expected_text)
            actual = extract_section(test_md, section_name)

            function project_for_equivalence(str)
                res = replace(str, r"(^|\n)\s*(\*|-)" => "*")
                res = strip(lowercase(replace(res, "\n" => " ")))
                while contains(res, "  ")
                    res = replace(res, "  " => " ")
                end
                res
            end

            expected_projected = project_for_equivalence(expected_text)
            actual_projected = project_for_equivalence(actual)

            function display_test_failure(expected, actual)
                println("The two strings, transformed for equality check:")
                println("Expected:")
                println(expected)
                println("Actual:")
                println(actual)
                println("Are not equal.")
            end

            if expected_projected != actual_projected
                display_test_failure(expected_projected, actual_projected)
            end

            @test expected_projected == actual_projected

        end

        map(
            pair -> test_section_extraction(pair.first, pair.second),
            [
                "Header 1" => "This is some text",
                "Header 2" => "This is some other text. On more paragraphs. This is a line in a paragraph.",
                "Header 3" => "This is the end bit. Which is quite long and has two paragraphs (and is at the end).",
                "Header 2.5" => """A list
                                   - one
                                   - two
                                   - three
                                """,
                "Header 2.6" => """A list
                                   * one
                                   * two
                                   * three
                                """,
                "A link" => "[this is a link](to somewhere far away)",
                "Some Code" => "```C int x,y ```",
            ],
        )

    end
end
