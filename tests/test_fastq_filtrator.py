import os.path

import pytest

from main import filter_fastq


def get_output_data(output):
    result = []
    with open(output) as file:
        for line in file:
            result.append(line.strip("\n"))
    return result


@pytest.fixture
def tmp_input():
    file_path = "tmp_input.fastq"
    with open(file_path, "w") as fastq:
        fastq.write("@one\nAAAAA\n+\n@@@@@\n")
        fastq.write("@two\nATGC\n+\n@@@@\n")
        fastq.write("@three\nGGG\n+\nEEE\n")

    yield file_path

    if os.path.exists(file_path):
        os.remove(file_path)


@pytest.fixture
def tmp_incorrect_fastq_input():
    file_path = "tmp_input_incorrect.fastq"
    with open(file_path, "w") as fastq:
        fastq.write("@incorr\n\n+\n\n")  # нет последоват-ти
        fastq.write("@corr\nGGG\n+\nFFF\n")

    yield file_path

    if os.path.exists(file_path):
        os.remove(file_path)


@pytest.fixture
def tmp_output():
    file_path = "tmp_output.fastq"

    yield file_path

    if os.path.exists(file_path):
        os.remove(file_path)


def test_filter_fastq_input_file_not_exists(tmp_output):
    with pytest.raises(FileNotFoundError):
        filter_fastq(input_fastq="....", output_fastq=tmp_output)


def test_filter_fastq_file_incorrect(tmp_incorrect_fastq_input, tmp_output):
    with pytest.raises(ValueError):
        filter_fastq(input_fastq=tmp_incorrect_fastq_input, output_fastq=tmp_output)


def test_filter_fastq_output_exists(tmp_input, tmp_output):
    filter_fastq(input_fastq=tmp_input, output_fastq=tmp_output)
    assert os.path.exists(tmp_output)


def test_filter_fastq_with_default(tmp_input, tmp_output):
    filter_fastq(input_fastq=tmp_input, output_fastq=tmp_output)

    expected_result = [
        "@one",
        "AAAAA",
        "+",
        "@@@@@",
        "@two",
        "ATGC",
        "+",
        "@@@@",
        "@three",
        "GGG",
        "+",
        "EEE",
    ]

    assert get_output_data(tmp_output) == expected_result


def test_filter_fastq_by_len(tmp_input, tmp_output):
    filter_fastq(input_fastq=tmp_input, output_fastq=tmp_output, length_bounds=5)

    expected_result = [
        "@one",
        "AAAAA",
        "+",
        "@@@@@",
    ]

    assert get_output_data(tmp_output) == expected_result


def test_filter_fastq_by_gc(tmp_input, tmp_output):
    filter_fastq(input_fastq=tmp_input, output_fastq=tmp_output, gc_bounds=(90, 100))

    expected_result = [
        "@three",
        "GGG",
        "+",
        "EEE",
    ]

    assert get_output_data(tmp_output) == expected_result


def test_filter_fastq_by_gc_only_lower_bound(tmp_input, tmp_output):
    filter_fastq(input_fastq=tmp_input, output_fastq=tmp_output, gc_bounds=90)

    expected_result = [
        "@three",
        "GGG",
        "+",
        "EEE",
    ]

    assert get_output_data(tmp_output) == expected_result


def test_filter_fastq_by_quality(tmp_input, tmp_output):
    filter_fastq(input_fastq=tmp_input, output_fastq=tmp_output, quality_threshold=32)

    expected_result = [
        "@three",
        "GGG",
        "+",
        "EEE",
    ]

    assert get_output_data(tmp_output) == expected_result
