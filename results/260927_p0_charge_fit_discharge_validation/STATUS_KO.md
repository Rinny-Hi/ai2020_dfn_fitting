# 중단 상태

사용자 요청에 따라 5개 파라미터 동시 fitting을 최종 결과 생성 전에 중단했다.

- `optimizer_evaluation_history.csv`는 중단 전 임시 평가 이력이며 최종 fitting
  결과로 사용하지 않는다.
- 이후 fitting은 mixed LSA/correlation을 통과한 `Dsn`, `kn`, `brugg_n` subset과
  비교용 축소 subset만 대상으로 새로 실행해야 한다.
- 이 폴더에는 확정된 최적 파라미터 또는 validation 결과가 없다.

