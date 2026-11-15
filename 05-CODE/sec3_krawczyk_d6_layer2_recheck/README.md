# README — sec3_krawczyk_d6_layer2_recheck

Пакет перепроверки двух боксов Кравчика слоя 2 (S078, S090) на семействе `1 − a·x⁶`, m = 9, период n = 1024. Источник чисел: `S090_S078_recert.json` и `REPORT-R2-S090.md`. Воспроизводимая среда: `prec_bits ∈ {256, 512, 1024}`, `radius ∈ {1e-30, 1e-31, 1e-32, 1e-33, 1e-34}`, два порядка вычисления степени (`pow`, `mul = x⁵·x`).

## Содержимое

| файл | назначение |
|---|---|
| `r2_sweep.py` | свип `(prec × radius × pow_order)` для S078, S090 + контролей; выдаёт таблицу `kappa`, `relwidth(F')`, `pass` |
| `r2_recert.py` | финальная пересертификация S078, S090 при `r = 1e-34`, `prec = 512`, `mul`; пишет `S090_S078_recert.json` |
| `S090_S078_recert.json` | recert-артефакт: sha256-пины, корни, метрики (`kappa`, `relwidth`, `headroom`), периоды, cross-variant pow/mul |
| `r2_table.json` | сводная таблица свипа; загружается как валидационный артефакт |
| `krawczyk6_certificate_local.py` | локальная копия генератора Krawczyk+Newton сертификата для семейства `1 − a·x⁶` |
| `check_pairwise_isolation.py` | контроль попарной изоляции боксов (защита от «склейки» соседних корней при уменьшении радиуса) |

## Воспроизведение

Из корня проекта (`wave5_stage3/`):

```bash
# 1. Свип prec × radius × pow_order (генерирует r2_table.json):
python 05-CODE/sec3_krawczyk_d6_layer2_recheck/r2_sweep.py \
    --radii 1e-30 1e-31 1e-32 1e-33 1e-34 \
    --precs 256 512 1024 \
    --orders pow mul \
    --out 05-CODE/sec3_krawczyk_d6_layer2_recheck/r2_table.json

# 2. Финальный recert (записывает S090_S078_recert.json):
python 05-CODE/sec3_krawczyk_d6_layer2_recheck/r2_recert.py \
    --radius 1e-34 \
    --prec 512 \
    --order mul \
    --out 05-CODE/sec3_krawczyk_d6_layer2_recheck/S090_S078_recert.json

# 3. Контроль попарной изоляции:
python 05-CODE/sec3_krawczyk_d6_layer2_recheck/check_pairwise_isolation.py
```

## Что показано

| режим | S078 (κ, relwidth\|F'\|, pass) | S090 (κ, relwidth\|F'\|, pass) |
|---|---|---|
| `r = 1e-30, prec = 256, pow` | pass, но relwidth \|F'\| ≈ 1.3 % | pass, relwidth \|F'\| ≈ 5.5 % |
| `r = 1e-31, prec = 256, pow` | pass | pass, relwidth ≈ 6.83 %, κ ≈ 0.09 |
| `r = 1e-31, prec = 256, mul` | pass | **fail**: 0 ∈ F'(X) |
| `r = 1e-32, prec = 256, mul` | pass, relwidth 15.5 % | pass, relwidth 15.5 % |
| `r = 1e-33, prec = 256, mul` | pass, relwidth 1.4 % | pass, relwidth 1.4 % |
| **`r = 1e-34, prec = 512, mul`** | **pass, relwidth 1.14e-7, κ 1.61e-7** | **pass, relwidth 1.36e-3, κ 1.93e-3** |

**Хрупкость при `r = 1e-31` в порядке `x⁵·x`**: для S090 шар F'(X) в mul-пути содержит нуль (строгая сертификация в этом порядке невозможна). Это воспроизведено независимой проверкой как 286/287 pass через mul-путь (см. `REPORT-R2-S090.md` §1).

**Комфортный запас при `r = 1e-34`**: один и тот же бокс проходит в ОБОИХ порядках (`pow` и `mul`) при `prec ∈ {256, 512, 1024}`. Критерии строгого контроля (i)–(iv) выполнены:
- (i) `K(X) ⊂ X` (headroom ≈ r для S078, ≈ 0.998·r для S090);
- (ii) `0 ∉ F'(X)` (lb ≈ 8.12·10¹⁴ для S078, ≈ 1.08·10¹⁸ для S090);
- (iii) периоды `inf |x₂ⱼ|` ≥ 0.053 на j = 0..9;
- (iv) пин к stored seed `|c − sto|_ub = 1.74e-24` (S078) / `1.85e-24` (S090), оба < 5e-24.

## Числа для B.2.2.1 (Appendix B)

В `03-APPENDIX-B-DATA-v06.md` §B.2.2.1 для S078/S090 использован именно радиус `r = 1e-34` (prec = 512, mul), а не авторский `r = 1e-31`. Альтернативный путь (наклонная арифметика / форма среднего значения для F') не понадобился: уменьшение радиуса при `prec = 512` уже даёт комфортный запас по всем четырём критериям.

## Воспроизводимость

- SHA256-пины зависимостей (`refined_seeds6.json`, `krawczyk6_certificates.json`, `calibration6.json`) записаны в `S090_S078_recert.json` поле `sha256_pins`; при их изменении recert должен быть повторён.
- Шапка отчёта `REPORT-R2-S090.md` содержит сводный итог свипа; полная таблица — `r2_table.json`.
- Все цифры выше взяты посимвольно из этих JSON.
