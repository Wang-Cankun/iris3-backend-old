import { UserExistGuard } from './user-exist.guard';

describe('UserExistGuard', () => {
  it('should be defined', () => {
    expect(new UserExistGuard()).toBeDefined();
  });
});
