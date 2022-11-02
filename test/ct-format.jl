using RNAstructure: parse_ct_format

ct_str1 =
"""
  5
    1 G       0    2    5    1
    2 A       1    3    0    2
    3 A       2    4    0    3
    4 A       3    5    0    4
    5 C       4    0    1    5
"""

ct_str2 =
""" 5   Title
    1 G       0    2    5    1
    2 A       1    3    0    2
    3 A       2    4    0    3
    4 A       3    5    0    4
    5 C       4    0    1    5
6 Another title
    1 C       0    2    5    1
    2 A       1    3    0    2
    3 A       2    4    0    3
    4 A       3    5    0    4
    5 G       4    6    1    5
    6 A       5    0    0    6
4 A pseudoknot structure
    1 G       0    2    3    1
    2 C       1    3    4    2
    3 C       2    4    1    3
    4 G       3    0    2    4
"""

ct_str3 = """
 1  A title with  two spaces
    1 G       0    0    0    1
"""

@testset "parse_ct_format" begin
    results = parse_ct_format(ct_str1)
    @test length(results) == 1
    title, seq, pt = results[1]
    @test title == ""
    @test seq == ["G", "A", "A", "A", "C"]
    @test pt == [5, 0, 0, 0, 1]

    results = parse_ct_format(ct_str2)
    @test length(results) == 3
    title, seq, pt = results[1]
    @test title == "Title"
    @test seq == ["G", "A", "A", "A", "C"]
    @test pt == [5, 0, 0, 0, 1]
    title, seq, pt = results[2]
    @test title == "Another title"
    @test seq == ["C", "A", "A", "A", "G", "A"]
    @test pt == [5, 0, 0, 0, 1, 0]
    title, seq, pt = results[3]
    @test title == "A pseudoknot structure"
    @test seq == ["G", "C", "C", "G"]
    @test pt == [3, 4, 1, 2]

    results = parse_ct_format(ct_str3)
    @test length(results) == 1
    title, seq, pt = results[1]
    @test title == "A title with  two spaces"
end
