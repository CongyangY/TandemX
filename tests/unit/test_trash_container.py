import pytest

from benchmarks.scripts.run_trash_container import cgroup_cpu, linux_time, native_products, tool_command


def test_container_resource_units_and_failure_status_are_preserved(tmp_path):
    path = tmp_path/'time.tsv'
    header = 'elapsed_seconds\tmax_rss_kib\tuser_seconds\tsystem_seconds\texit_code\n'
    path.write_text(header+'Command exited with non-zero status 7\n1.25\t2048\t0.75\t0.1\t7\n')
    measured = linux_time(path)
    assert measured['peak_rss_mib'] == 2
    assert measured['exit_code'] == 7
    assert measured['elapsed_seconds'] == 1.25
    for body in ('', 'nan\t2048\t0\t0\t0\n', '1\t2048\t0\t0\t0.5\n',
                 '1\t2048\t0\t0\t0\n1\t2048\t0\t0\t0\n'):
        path.write_text(header+body)
        with pytest.raises(ValueError):
            linux_time(path)


def test_native_comparator_commands_use_distinct_author_interfaces_without_templates():
    assert tool_command('trash')[-4:] == ['--par', '1', '--randomseed', '6101']
    assert tool_command('trash2')[-2:] == ['-p', '1']
    assert '--def' in tool_command('trash')
    assert all('template' not in arg for t in ('trash', 'trash2') for arg in tool_command(t))
    assert native_products('trash2')[1] == ('assembly.fa_arrays.csv', 'assembly.fa_repeats.csv')
    with pytest.raises(ValueError):
        tool_command('other')


def test_cgroup_cpu_includes_all_tasks_and_preserves_counter_units(tmp_path):
    before, after = tmp_path/'before', tmp_path/'after'
    before.write_text('usage_usec 100\nuser_usec 80\nsystem_usec 20\nnr_throttled 2\n')
    after.write_text('usage_usec 2500100\nuser_usec 2000080\nsystem_usec 500020\nnr_throttled 12\n')
    measured = cgroup_cpu(before, after)
    assert measured['usage_seconds'] == 2.5
    assert measured['user_seconds'] == 2
    assert measured['system_seconds'] == .5
    assert measured['delta_counters']['nr_throttled'] == 10
    with pytest.raises(ValueError, match='decreased'):
        cgroup_cpu(after, before)
    after.write_text('usage_usec 100\nusage_usec 100\n')
    with pytest.raises(ValueError, match='Invalid'):
        cgroup_cpu(before, after)
