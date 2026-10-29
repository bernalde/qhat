"""User-facing model parsing and certificate disclosure, without simulations."""
import json
import pytest

from qhat.analysis.train_ordering import load_model, main


def model_file(tmp_path):
    path = tmp_path / 'model.json'
    path.write_text(json.dumps(dict(n=1, electrons=0, denominator='2',
                                    blocks=[[[[[0, 'X']], '2']], [[[[0, 'Z']], '1']]])))
    return path


def test_command_reports_real_order_and_separate_certificate(tmp_path, capsys):
    path = model_file(tmp_path)
    assert main([str(path), '--rotation-budget', '60', '--cpu-seconds', '1',
                 '--max-evaluations', '3', '--error-tolerance', '.01']) == 0
    data = json.loads(capsys.readouterr().out)
    assert sorted(data['training']['order']) == [0, 1]
    assert data['training']['score']['value'] >= 0
    assert data['certificate']['status'] == 'certified_ideal_bound'
    assert data['certificate']['order_dependent'] is False
    assert 'not a finite-error bound' in data['score_scope']
    assert 'hardware' in data['certificate_scope']


def test_command_discloses_unresolved_step_cap_and_refuses_overwrite(tmp_path):
    path = model_file(tmp_path)
    output = tmp_path / 'result.json'
    args = [str(path), '--rotation-budget', '30', '--max-evaluations', '0',
            '--error-tolerance', '1e-8', '--max-certified-steps', '1', '--output', str(output)]
    assert main(args) == 0
    data = json.loads(output.read_text())
    assert data['training']['score'] is None
    assert data['certificate']['status'] == 'unresolved_within_step_limit'
    with pytest.raises(FileExistsError):
        main(args)


@pytest.mark.parametrize('field,bad', [('denominator', True), ('denominator', 1.5)])
def test_denominator_is_not_silently_coerced(tmp_path, field, bad):
    path = model_file(tmp_path)
    model = json.loads(path.read_text());model[field] = bad;path.write_text(json.dumps(model))
    with pytest.raises(ValueError, match='integers'):
        load_model(path)


def test_coefficient_is_not_silently_truncated(tmp_path):
    path = model_file(tmp_path)
    model = json.loads(path.read_text());model['blocks'][0][0][1] = 1.5;path.write_text(json.dumps(model))
    with pytest.raises(ValueError, match='integers'):
        load_model(path)
