"""One exact weighted-quadratic coordinate pass for diagonal-term placement."""
from __future__ import annotations


def score(amplitudes, weights):
    return sum(weights[d]*a*a for d,a in amplitudes.items() if a)


def move(amplitudes, contributions, masks, old_slot, weights):
    trial=dict(amplitudes)
    for j in range(old_slot):
        d=masks[j];trial[d]=trial.get(d,0)+2*contributions[j]
    value=score(trial,weights)
    values=[value]
    for d,b in zip(masks,contributions):
        if b:
            before=trial.get(d,0);after=before-2*b
            value+=weights[d]*(after*after-before*before)
            trial[d]=after
        values.append(value)
    selected=min(range(len(values)),key=lambda slot:values[slot])
    result=dict(amplitudes)
    for j,(d,b) in enumerate(zip(masks,contributions)):
        change=((1 if j>=selected else -1)-(1 if j>=old_slot else -1))*b
        result[d]=result.get(d,0)+change
    assert score(result,weights)==values[selected]
    assert values[selected]<=score(amplitudes,weights)
    return selected,result,values


def place(off_amplitudes, contributions, masks, start_slot, weights):
    """Rows of contributions must already follow coefficient/word priority."""
    slots=[start_slot]*len(contributions)
    amplitudes=dict(off_amplitudes)
    for row in contributions:
        for j,(d,b) in enumerate(zip(masks,row)):
            amplitudes[d]=amplitudes.get(d,0)+(b if j>=start_slot else -b)
    history=[]
    for index,row in enumerate(contributions):
        before=score(amplitudes,weights)
        slot,amplitudes,values=move(amplitudes,row,masks,slots[index],weights)
        history.append(dict(index=index,old_slot=slots[index],new_slot=slot,
                            old_score=str(before),new_score=str(values[slot]),
                            slot_scores=[str(v) for v in values]))
        slots[index]=slot
    return slots,amplitudes,history
